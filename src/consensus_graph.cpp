#include <deps/odgi/src/odgi.hpp>
#include "consensus_graph.hpp"

namespace smoothxg {

bool operator<(const link_path_t& a,
               const link_path_t& b) {
    auto& a_0 = as_integer(a.from_cons_path);
    auto& a_1 = as_integer(a.to_cons_path);
    auto& b_0 = as_integer(b.from_cons_path);
    auto& b_1 = as_integer(b.to_cons_path);
    return (a_0 < b_0) 
        || (a_0 == b_0
            && (a_1 < b_1
                || (a_1 == b_1
                    && (a.length < b.length
                        || (a.length == b.length
                            && a.hash < b.hash)))));
}

ostream& operator<<(ostream& o, const link_path_t& a) {
    o << "("
      << a.from_cons_name << " "
      << a.to_cons_name << " "
      << as_integer(a.from_cons_path) << " "
      << as_integer(a.to_cons_path) << " "
      << a.length << " "
      << a.hash << " "
      << as_integer(a.path) << " "
      << as_integers(a.begin)[0] << ":" << as_integers(a.begin)[1] << " "
      << as_integers(a.end)[0] << ":" << as_integers(a.end)[1] << ")";
    return o;
}

std::vector<consensus_spec_t> parse_consensus_spec(const std::string& spec_str,
                                                   bool& requires_consensus) {
    auto fields = split(spec_str, ',');
    assert(fields.size() > 1);
    auto& basename = fields.front();
    std::vector<consensus_spec_t> specs;
    for (int i = 1; i < fields.size(); ++i) {
        // split the field, add the spec
        specs.emplace_back();
        auto& spec = specs.back();
        spec.basename = basename;
        auto vals = split(fields[i], ':');
        if (vals.size() > 0) {
            spec.min_allele_len = std::stoi(vals[0]);
        }
        if (vals.size() > 1) {
            spec.ref_file = vals[1];
            // sanitize file name
            spec.ref_file_sanitized = spec.ref_file;
            for (auto& c : spec.ref_file_sanitized) {
                if (c == '/') c = '_';
            }
        }
        if (vals.size() > 2) {
            spec.keep_consensus_paths = (vals[2] == "y");
        } else {
            spec.keep_consensus_paths = true;
        }
        if (vals.size() > 3) {
            spec.min_consensus_path_cov = std::stod(vals[3]);
        } else {
            spec.min_consensus_path_cov = 0;
        }
        if (vals.size() > 4) {
            spec.max_allele_len = std::stoi(vals[4]);
        } else {
            spec.max_allele_len = 1e6;
        }

        requires_consensus |= spec.keep_consensus_paths;
    }
    return specs;
}

std::string displayname(const consensus_spec_t& spec) {
    std::stringstream s;
    s << spec.basename << "@"
      << spec.min_allele_len
      << "_" << (!spec.ref_file.empty() ? spec.ref_file_sanitized : "")
      << "_" << (spec.keep_consensus_paths ? "y" : "n")
      << "_" << spec.min_consensus_path_cov
      << "_" << spec.max_allele_len;
    return s.str();
}


// prep the graph into a given GFA file
// we'll then build the xg index on top of that in low memory

odgi::graph_t* create_consensus_graph(const xg::XG &smoothed,
                                      // TODO: GBWT
                                      const std::vector<std::string>& _consensus_path_names,
                                      const uint64_t& min_allele_length,
                                      const uint64_t& max_allele_length,
                                      const double& min_consensus_path_coverage,
                                      // TODO: minimum allele frequency
                                      const uint64_t& thread_count,
                                      const std::string& base) {



    // OVERALL: https://www.acodersjourney.com/6-tips-supercharge-cpp-11-vector-performance/ -> check these things here

    std::vector<path_handle_t> consensus_paths;
    //consensus_paths.reserve(consensus_path_names.size());
    std::vector<std::string> consensus_path_names;
    for (auto& name : _consensus_path_names) {
        if (smoothed.has_path(name)) {
            consensus_paths.push_back(smoothed.get_path_handle(name));
            consensus_path_names.push_back(name);
        }
    }
    if (consensus_paths.empty()) {
        std::cerr << "[smoothxg::create_consensus_graph] WARNING: no matching paths found, returning an empty graph" << std::endl;
        return new odgi::graph_t();
    }

    // compute mean coverage for the paths if we are filtering
    // we will remove consensus paths lower than this coverage
    if (min_consensus_path_coverage) {
        std::vector<double> consensus_mean_coverage(consensus_paths.size());
#pragma omp parallel for schedule(static, 1) num_threads(thread_count)
        for (uint64_t i = 0; i < consensus_paths.size(); ++i) {
            uint64_t length = 0;
            uint64_t coverage = 0;
            smoothed.for_each_step_in_path(
                consensus_paths[i],
                [&](const step_handle_t& step) {
                    handle_t handle = smoothed.get_handle_of_step(step);
                    uint64_t handle_length = smoothed.get_length(handle);
                    length += handle_length;
                    uint64_t depth = 0;
                    smoothed.for_each_step_on_handle(
                        handle, [&](const step_handle_t& s) { ++depth; });
                    coverage += length * depth;
                });
            consensus_mean_coverage[i] = (double)coverage / (double)length;
        }
        std::vector<path_handle_t> cov_consensus_paths;
        std::vector<std::string> cov_consensus_path_names;
        for (uint64_t i = 0; i < consensus_paths.size(); ++i) {
            if (consensus_mean_coverage[i] > min_consensus_path_coverage) {
                cov_consensus_paths.push_back(consensus_paths[i]);
                cov_consensus_path_names.push_back(smoothed.get_path_name(consensus_paths[i]));
            }
        }
        consensus_paths = cov_consensus_paths;
        consensus_path_names = cov_consensus_path_names;
    }

    std::vector<bool> is_consensus(smoothed.get_path_count()+1, false);
    for (auto& path : consensus_paths) {
        is_consensus[as_integer(path)] = true;
    }

    std::vector<std::string*> cons_path_ptr(smoothed.get_path_count()+1, nullptr);
    uint64_t idx = 0;
    for (auto& path : consensus_paths) {
        cons_path_ptr[as_integer(path)] = (std::string*)&consensus_path_names[idx++];
    }

    // here we compute a consensus path per node
    // TODO we should allow a vector or set of consensus paths per node
    // and we mark a vector that says if we're in a consensus or not
    std::vector<path_handle_t> consensus_path_handles(smoothed.get_node_count());
    atomicbitvector::atomic_bv_t handle_is_consensus(smoothed.get_node_count());
#pragma omp parallel for schedule(static, 1) num_threads(thread_count)
    for (uint64_t i = 0; i < consensus_paths.size(); ++i) {
        uint64_t length = 0;
        uint64_t coverage = 0;
        smoothed.for_each_step_in_path(
            consensus_paths[i],
            [&](const step_handle_t& step) {
                handle_t handle = smoothed.get_handle_of_step(step);
                nid_t node_id = smoothed.get_id(handle);
                // save a consensus path for each, the first we get to
                if (!handle_is_consensus.set(node_id - 1)) {
                    consensus_path_handles[node_id - 1] = consensus_paths[i];
                }
            });
    }

    std::vector<path_handle_t> non_consensus_paths;
    non_consensus_paths.reserve(smoothed.get_path_count()+1-consensus_paths.size());
    smoothed.for_each_path_handle(
        [&](const path_handle_t& p) {
            if (!is_consensus[as_integer(p)]) {
                non_consensus_paths.push_back(p);
            }
        });

    auto get_path_seq_length =
        [&](const step_handle_t& begin,
            const step_handle_t& end) {
            uint64_t len = 0;
            for (step_handle_t i = begin; i != end; i = smoothed.get_next_step(i)) {
                len += smoothed.get_length(smoothed.get_handle_of_step(i));
            }
            return len;
        };
    
    auto get_path_seq =
        [&](const step_handle_t& begin,
            const step_handle_t& end) {
            std::string seq;
            for (step_handle_t i = begin; i != end; i = smoothed.get_next_step(i)) {
                seq.append(smoothed.get_sequence(smoothed.get_handle_of_step(i)));
            }
            return seq;
        };

    auto hash_seq =
        [&](const std::string& seq) {
            return std::hash<std::string>{}(seq);
        };

    auto start_in_vector =
        [&](const handle_t& h) {
            if (!smoothed.get_is_reverse(h)) {
                return (int64_t) smoothed.node_vector_offset(smoothed.get_id(h));
            } else {
                return (int64_t) (smoothed.node_vector_offset(smoothed.get_id(h))
                                  + smoothed.get_length(h));
            }
        };

    auto end_in_vector =
        [&](const handle_t& h) {
            if (smoothed.get_is_reverse(h)) {
                return (int64_t) smoothed.node_vector_offset(smoothed.get_id(h));
            } else {
                return (int64_t) (smoothed.node_vector_offset(smoothed.get_id(h))
                                  + smoothed.get_length(h));
            }
        };

    // consensus path -> consensus path : link_path_t
    std::string base_mmset = base + ".link_path_ms";
    auto link_path_ms = std::make_unique<mmmulti::set<link_path_t>>(base_mmset);
    link_path_ms->open_writer();

    // TODO: parallelize over path ranges that tend to have around the same max length
    // determine the ranges based on a map of the consensus path set
    std::atomic<bool> is_there_something(false);
    // TODO: this could reflect the haplotype frequencies to preserve variation > some frequency
#pragma omp parallel for schedule(static, 1) num_threads(thread_count)
    for (uint64_t idx = 0; idx < non_consensus_paths.size(); ++idx){
        auto& path = non_consensus_paths[idx];
        // for each step in path
        link_path_t link;
        link.path = path;
        path_handle_t last_seen_consensus;
        bool seen_consensus = false;
        smoothed.for_each_step_in_path(path, [&](const step_handle_t& step) {
            // check if we're on the step with any consensus
            handle_t h = smoothed.get_handle_of_step(step);
            nid_t node_id = smoothed.get_id(h);
            bool on_consensus = false;
            path_handle_t curr_consensus;
            // we use a bitvector (vector<bool>) saying if the node is in the consensus graph
            // and then we don't have to iterate path_steps * path_steps
            // we'll also need to store the path handle of the consensus path at that node (another vector!)
            // .... but keep in mind that this makes the assumption that we have only one consensus path at any place in the graph
            // .... if we have more, should we just handle the first in this context...?
            ///.... currently we are using the "last" one we find (but there's only one)
            if (handle_is_consensus.test(node_id - 1)) {
                on_consensus = true;
                curr_consensus = consensus_path_handles[node_id - 1];
            }
            // if we're on the consensus
            if (on_consensus) {
                // we haven't seen any consensus before?
                if (!seen_consensus) {
                    // we construct the first link path object
                    link.length = 0;
                    link.from_cons_name = cons_path_ptr[as_integer(curr_consensus)];
                    link.to_cons_name = cons_path_ptr[as_integer(curr_consensus)];
                    link.from_cons_path = curr_consensus;
                    link.to_cons_path = curr_consensus;
                    link.begin = step;
                    link.end = step;
                    link.hash = 0;
                    seen_consensus = true;
                    last_seen_consensus = curr_consensus;
                    // TODO do we want to add the allele depth of the start and end consensus handles to the link object?

                    // TODO add frequency to construct link_path_t --> double
                    // TODO min and max frequency
                } else {
                    /*
                    if (last_seen_consensus != consensus) {
                        std::cerr << "path " << smoothed.get_path_name(path) << " switched from " << smoothed.get_path_name(last_seen_consensus) << " to " << smoothed.get_path_name(consensus) << std::endl;
                        last_seen_consensus = consensus;
                    }
                    */

                    // we've seen a consensus before, and it's the same
                    // and the direction of movement is correct
                    // check the distance in the graph position vector
                    // if it's over some threshold, record the link
                    handle_t last_handle = smoothed.get_handle_of_step(link.end);
                    handle_t curr_handle = smoothed.get_handle_of_step(step);
                    uint64_t jump_length = std::abs(start_in_vector(curr_handle)
                                                    - end_in_vector(last_handle));

                    // TODO: don't just look at min_allele_length, but also consider allele frequency
                    if (link.from_cons_path == curr_consensus
                        && jump_length < min_allele_length) {
                        link.begin = step;
                        link.end = step;
                        link.length = 0;
                    } else { // or it's different
                        // this is when we write a link candidate record
                        link.to_cons_name = cons_path_ptr[as_integer(curr_consensus)];
                        link.to_cons_path = curr_consensus;
                        //link.begin = smoothed.get_next_step(link.begin);
                        link.end = step;
                        //std::cerr << "writing to mmset" << std::endl;
                        link.length = get_path_seq_length(
                                smoothed.get_next_step(link.begin),
                                link.end);
                        stringstream h;
                        std::string seq = get_path_seq(smoothed.get_next_step(link.begin),
                                                       link.end);
                        h << smoothed.get_id(smoothed.get_handle_of_step(link.begin))
                          << ":"
                          << smoothed.get_id(smoothed.get_handle_of_step(link.end))
                          << ":"
                          << seq;
                        link.hash = hash_seq(h.str());
                        if (as_integer(link.from_cons_path) > as_integer(link.to_cons_path)) {
                            std::swap(link.from_cons_path, link.to_cons_path);
                            std::swap(link.from_cons_name, link.to_cons_name);
                        }
                        link.jump_length = jump_length;
                        link_path_ms->append(link);
                        is_there_something.store(true);

                        // reset link
                        link.length = 0;
                        link.from_cons_name = cons_path_ptr[as_integer(curr_consensus)];
                        link.to_cons_name = cons_path_ptr[as_integer(curr_consensus)];
                        link.from_cons_path = curr_consensus;
                        link.to_cons_path = curr_consensus;
                        link.begin = step;
                        link.end = step;
                        link.hash = 0;
                        link.is_rev = smoothed.get_is_reverse(smoothed.get_handle_of_step(step));
                    }
                }
            } else {
                //link.length += smoothed.get_length(h);
            }
        });
    }

    std::vector<link_path_t> consensus_links;
    std::vector<std::pair<handle_t, handle_t>> perfect_edges;

    if (is_there_something.load()) {
        link_path_ms->index(thread_count);

        // collect sets of link paths that refer to the same consensus path pairs
        // and pick which one to keep in the consensus graph

        std::vector<link_path_t> curr_links;

        path_handle_t curr_from_cons_path;
        path_handle_t curr_to_cons_path;

        auto novel_sequence_length =
                [&](const step_handle_t begin,
                    const step_handle_t end,
                    ska::flat_hash_set<uint64_t>& seen_nodes, // by ref
                    const xg::XG &graph) {
                    uint64_t novel_bp = 0;
                    for (auto s = begin;
                         s != end; s = graph.get_next_step(s)) {
                        handle_t h = graph.get_handle_of_step(s);
                        uint64_t i = graph.get_id(h);
                        if (!seen_nodes.count(i)) {
                            novel_bp += graph.get_length(h);
                            seen_nodes.insert(i);
                        }
                    }
                    return novel_bp;
                };

        auto largest_novel_gap =
                [&](const step_handle_t begin,
                    const step_handle_t end,
                    ska::flat_hash_set<uint64_t>& seen_nodes, // by ref
                    const xg::XG &graph) {
                    uint64_t novel_bp = 0;
                    uint64_t largest_gap = 0;
                    for (auto s = begin;
                         s != end; s = graph.get_next_step(s)) {
                        handle_t h = graph.get_handle_of_step(s);
                        uint64_t i = graph.get_id(h);
                        if (!seen_nodes.count(i)) {
                            novel_bp += graph.get_length(h);
                            seen_nodes.insert(i);
                        } else {
                            largest_gap = std::max(novel_bp, largest_gap);
                            novel_bp = 0;
                        }
                    }
                    return largest_gap;
                };

        auto get_step_count =
                [&](const step_handle_t begin,
                    const step_handle_t end,
                    const xg::XG &graph) {
                    uint64_t count = 0;
                    for (auto s = begin;
                         s != end; s = graph.get_next_step(s)) {
                        ++count;
                    }
                    return count;
                };


        auto mark_seen_nodes =
                [&](const step_handle_t begin,
                    const step_handle_t end,
                    ska::flat_hash_set<uint64_t> &seen_nodes, // by ref
                    const xg::XG &graph) {
                    for (auto s = begin;
                         s != end; s = graph.get_next_step(s)) {
                        handle_t h = graph.get_handle_of_step(s);
                        uint64_t i = graph.get_id(h);
                        if (!seen_nodes.count(i)) {
                            seen_nodes.insert(i);
                        }
                    }
                };

        // FIXME can we parallelize this?
        auto compute_best_link =
                [&](const std::vector<link_path_t> &links) {
                    std::map<uint64_t, uint64_t> hash_counts;
                    std::vector<link_path_t> unique_links;
                    uint64_t link_rank = 0;
                    for (auto &link : links) {
                        //std::cerr << link << std::endl;
                        auto &c = hash_counts[link.hash];
                        if (c == 0) {
                            unique_links.push_back(link);
                        }
                        ++c;
                    }
                    std::map<uint64_t, uint64_t> hash_lengths;
                    for (auto &link : links) {
                        hash_lengths[link.hash] = link.length;
                    }
                    uint64_t best_count = 0;
                    uint64_t best_hash;
                    for (auto &c : hash_counts) {
                        if (c.second > best_count) {
                            best_hash = c.first;
                            best_count = c.second;
                        }
                    }
                    //std::cerr << "best hash be " << best_hash << std::endl;
                    // save the best link path
                    link_path_t most_frequent_link;
                    for (auto &link : unique_links) {
                        if (link.hash == best_hash) {
                            most_frequent_link = link;
                            break;
                        }
                    }

                    // if we have a 0-length link between consensus ends, add it
                    handle_t from_end_fwd
                            = smoothed.get_handle_of_step(
                                    smoothed.path_back(
                                            most_frequent_link.from_cons_path));
                    handle_t to_begin_fwd
                            = smoothed.get_handle_of_step(
                                    smoothed.path_begin(
                                            most_frequent_link.to_cons_path));
                    handle_t from_end_rev = smoothed.flip(to_begin_fwd);
                    handle_t to_begin_rev = smoothed.flip(from_end_fwd);

                    handle_t from_begin_fwd
                            = smoothed.get_handle_of_step(
                                    smoothed.path_begin(
                                            most_frequent_link.from_cons_path));
                    handle_t to_end_fwd
                            = smoothed.get_handle_of_step(
                                    smoothed.path_back(
                                            most_frequent_link.to_cons_path));
                    handle_t from_begin_rev = smoothed.flip(from_begin_fwd);
                    handle_t to_end_rev = smoothed.flip(to_end_fwd);

                    // if we can walk forward on the last handle of consensus a and reach consensus b
                    // we should add a link and say we're perfect

                    bool has_perfect_edge = false;
                    bool has_perfect_link = false;
                    link_path_t perfect_link;

                    if (smoothed.has_edge(from_end_fwd, to_begin_fwd)) {
                        auto p = std::make_pair(from_end_fwd, to_begin_fwd);
                        perfect_edges.push_back(p);
                        has_perfect_edge = true;
                    } else if (smoothed.has_edge(to_end_fwd, from_begin_fwd)) {
                        auto p = std::make_pair(to_end_fwd, from_begin_fwd);
                        perfect_edges.push_back(p);
                        has_perfect_edge = true;
                    } else {
                        for (auto &link : unique_links) {
                            //path_handle_t from_cons;
                            //path_handle_t to_cons;
                            // are we just stepping from the end of one to the beginning of the other?
                            for (step_handle_t s = link.begin;
                                 s != link.end; s = smoothed.get_next_step(s)) {
                                step_handle_t next = smoothed.get_next_step(s);
                                handle_t b = smoothed.get_handle_of_step(s);
                                handle_t e = smoothed.get_handle_of_step(next);
                                if (b == from_end_fwd && e == to_begin_fwd
                                    || b == from_end_rev && e == to_begin_rev
                                    || b == to_begin_fwd && e == from_end_fwd
                                    || b == to_begin_rev && e == from_end_rev) {
                                    has_perfect_link = true;
                                    perfect_link = link;
                                    break;
                                }
                            }
                            if (has_perfect_link) {
                                break;
                            }
                        }
                    }

                    ska::flat_hash_set<uint64_t> seen_nodes;

                    // this part attempts to preserve connectivity between consensus sequences
                    // we're preserving the consensus graph topology
                    if (has_perfect_edge) {
                        // nothing to do
                    } else if (has_perfect_link) {
                        // nb: should be no nodes to mark
                        mark_seen_nodes(perfect_link.begin, perfect_link.end, seen_nodes, smoothed);
                        perfect_link.rank = link_rank++;
                        consensus_links.push_back(perfect_link);
                    } else {
                        if (most_frequent_link.from_cons_path != most_frequent_link.to_cons_path) {
                            most_frequent_link.rank = link_rank++;
                            consensus_links.push_back(most_frequent_link);
                            mark_seen_nodes(most_frequent_link.begin, most_frequent_link.end, seen_nodes, smoothed);
                        }
                    }

                    // this part collects sequences that diverge from a consensus for the
                    // min_allele_length which is the "variant scale factor" of the algorithm
                    // this preserves novel non-consensus sequences greater than this length
                    // TODO: this could be made to respect allele frequency
                    for (auto &link : unique_links) {
                        if (link.hash == best_hash) {
                            continue;
                        }
                        // use the largest gap to determine if the link path is over our scale threshold
                        uint64_t largest_novel_gap_bp = largest_novel_gap(link.begin, link.end, seen_nodes, smoothed);
                        //uint64_t novel_bp = novel_sequence_length(link.begin, link.end, seen_nodes, smoothed);
                        uint64_t step_count = get_step_count(link.begin, link.end, smoothed);
                        // this complex filter attempts to keep representative link paths for indels above our min_allele_length
                        // we either need the jump length (measured in terms of delta in our graph vector) to be over our jump max
                        // *and* the link path should be empty or mostly novel
                        // *or* we're adding in the specified amount of novel_bp of sequence
                        if ((link.jump_length >= min_allele_length
                             && link.jump_length < max_allele_length
                             && (link.length == 0 ||
                                 (double)largest_novel_gap_bp / (double)step_count > 1))
                            || (largest_novel_gap_bp >= min_allele_length
                                && largest_novel_gap_bp < max_allele_length)) {
                            link.rank = link_rank++;
                            consensus_links.push_back(link);
                            mark_seen_nodes(link.begin, link.end, seen_nodes, smoothed);
                        }
                    }
                };

        std::cerr << "[smoothxg::create_consensus_graph] finding consensus links" << std::endl;
        // collect edges by node
        //
        // could this run in parallel?
        // FIXME this is already paralellized now, right?
        // yes probably but we need to either lock the consensus path vectors or write into a multiset
        link_path_ms->for_each_value(
                [&](const link_path_t &v) {
                    //std::cerr << "on " << v << " with count " << c << std::endl;
                    if (curr_links.empty()) {
                        curr_from_cons_path = v.from_cons_path;
                        curr_to_cons_path = v.to_cons_path;
                        curr_links.push_back(v);
                    } else if (curr_from_cons_path != v.from_cons_path
                               || curr_to_cons_path != v.to_cons_path) {
                        // compute the best link in the set of curr_links
                        compute_best_link(curr_links);
                        // reset the links
                        curr_links.clear();
                        curr_from_cons_path = v.from_cons_path;
                        curr_to_cons_path = v.to_cons_path;
                        curr_links.push_back(v);
                    } else {
                        curr_links.push_back(v);
                    }
                });

        compute_best_link(curr_links);
    }

    link_path_ms->close_reader();
    link_path_ms.reset();
    std::remove(base_mmset.c_str());

    // create new consensus graph which only has the consensus and link paths in it
    auto* consensus = new odgi::graph_t();
    consensus->set_number_of_threads(thread_count);

    // add the consensus paths first
    // FIXME could this be run in parallel? --> create the path in an extra for
    std::cerr << "[smoothxg::create_consensus_graph] adding consensus paths" << std::endl;
    for (auto& path : consensus_paths) {
        // create the path
        path_handle_t path_cons_graph = consensus->create_path_handle(smoothed.get_path_name(path));
        handle_t cur_handle_in_cons_graph;
        // add the current node first, then add the step
        smoothed.for_each_step_in_path(path,
                                       [&consensus, &smoothed, &path, &cur_handle_in_cons_graph, &path_cons_graph]
                                       (const step_handle_t& step) {
            handle_t h = smoothed.get_handle_of_step(step);
            bool rev = smoothed.get_is_reverse(h);
            nid_t node_id = smoothed.get_id(h);
            // we make a new node with the same sequence and relative orientation
            if (!consensus->has_node(node_id)) {
                if (!rev) {
                    cur_handle_in_cons_graph = consensus->create_handle(smoothed.get_sequence(h), node_id);
                } else {
                    cur_handle_in_cons_graph = consensus->create_handle(smoothed.get_sequence(smoothed.flip(h)), node_id);
                }
            } else {
                cur_handle_in_cons_graph = consensus->get_handle(node_id);
            }
            // our cur_handle_in_cons_graph is in forward orientation, so we have to match
            if (rev) {
                consensus->append_step(path_cons_graph, consensus->flip(cur_handle_in_cons_graph));
            } else {
                consensus->append_step(path_cons_graph, cur_handle_in_cons_graph);
            };
        });
    }

    std::cerr << "[smoothxg::create_consensus_graph] adding link paths" << std::endl;

    /// FIXME could THIS run in parallel
    // add link paths and edges not in the consensus paths
    std::vector<std::string> link_path_names;
    for (auto& link : consensus_links) {
        //std::cerr << "making " << "Link_" << smoothed.get_path_name(link.from_cons) << "_" << smoothed.get_path_name(link.to_cons) << "_" << link.rank << std::endl;
        // create link paths and paths
        // make the path name
        if (link.length > 0) {
            auto& novel_link = link;
            stringstream s;
            s << "Link_" << *novel_link.from_cons_name << "_" << *novel_link.to_cons_name << "_" << novel_link.rank;
            assert(!consensus->has_path(s.str()));
            path_handle_t path_cons_graph = consensus->create_path_handle(s.str());
            link_path_names.push_back(s.str());
            //link_paths.push_back(path_cons_graph);
            handle_t cur_handle_in_cons_graph;
            // add the current node first, then add the step
            for (step_handle_t step = novel_link.begin; //smoothed.get_next_step(novel_link.begin);
                 step != novel_link.end;
                 step = smoothed.get_next_step(step)) {
                handle_t h = smoothed.get_handle_of_step(step);
                handle_t next_handle;
                nid_t node_id = smoothed.get_id(h);
                if (!consensus->has_node(node_id)) {
                    //assert(false);
                    // this must be kept
                    cur_handle_in_cons_graph = consensus->create_handle(smoothed.get_sequence(smoothed.get_handle(node_id)), node_id);
                } else {
                    cur_handle_in_cons_graph = consensus->get_handle(node_id);
                }
                bool rev = smoothed.get_is_reverse(h);
                if (rev) {
                    consensus->append_step(path_cons_graph, consensus->flip(cur_handle_in_cons_graph));
                } else {
                    consensus->append_step(path_cons_graph, cur_handle_in_cons_graph);
                }
            }
        }
    }

    // finally add the edges
    consensus->for_each_path_handle(
        [&](const path_handle_t& path) {
            consensus->for_each_step_in_path(path, [&] (const step_handle_t step) {
               if (consensus->has_next_step(step)) {
                   step_handle_t next_step = consensus->get_next_step(step);
                   handle_t h = consensus->get_handle_of_step(step);
                   handle_t next_h = consensus->get_handle_of_step(next_step);
                   if (!consensus->has_edge(h, next_h)) {
                       consensus->create_edge(h, next_h);
                   }
               }
            });
        });

    for (auto& e : perfect_edges) {
        handle_t h = consensus->get_handle(
            smoothed.get_id(e.first),
            smoothed.get_is_reverse(e.first));
        handle_t j = consensus->get_handle(
            smoothed.get_id(e.second),
            smoothed.get_is_reverse(e.second));
        consensus->create_edge(h, j);
    }

    auto link_steps =
        [&](const step_handle_t& a, const step_handle_t& b) {
            handle_t from = smoothed.get_handle_of_step(a);
            handle_t to = smoothed.get_handle_of_step(b);
            nid_t from_id = smoothed.get_id(from);
            nid_t to_id = smoothed.get_id(to);
            if (consensus->has_node(from_id)
                && consensus->has_node(to_id)) {
                consensus->create_edge(consensus->get_handle(smoothed.get_id(from),
                                                           smoothed.get_is_reverse(from)),
                                      consensus->get_handle(smoothed.get_id(to),
                                                           smoothed.get_is_reverse(to)));
            }
        };

    // preserve topology of links by walking the path they derive from
    // forward and backward and ensuring this is in the consensus graph
    for (auto& link : consensus_links) {
        // edge from begin to begin+1
        // edge from end-1 to end
        step_handle_t next = smoothed.get_next_step(link.begin);
        link_steps(link.begin, next);
        step_handle_t prev = smoothed.get_previous_step(link.end);
        if (prev != link.begin) {
            link_steps(prev, link.end);
        }
    }

    /// TODO validate consensus graph until here
    // FIXME: this should check the actual path sequence for validation
    // not each step
    /*
    smoothed.for_each_path_handle(
            [&](const path_handle_t &p) {
                if (is_consensus[as_integer(p)]) {
                    std::string path_name = smoothed.get_path_name(p);
                    if (!consensus->has_path(path_name)) {
                        std::cerr << "[smoothxg::main::create_consensus_graph] error: consensus path " << path_name
                        << " not present in the consensus graph!" << std::endl;
                        exit(1);
                    }
                    smoothed.for_each_step_in_path(p,
                                                   [&]
                                                   (const step_handle_t& step) {
                        handle_t h = smoothed.get_handle_of_step(step);
                        nid_t node_id = smoothed.get_id(h);
                        // node id comparison
                        if (!consensus->has_node(node_id)) {
                            std::cerr << "[smoothxg::main::create_consensus_graph] error: node " << node_id
                            << " not present in the consensus graph!" << std::endl;
                            exit(1);
                        }
                        // node orientation comparison
                        handle_t consensus_h = consensus->get_handle(node_id, smoothed.get_is_reverse(h));
                        if (consensus->get_is_reverse(consensus_h) != smoothed.get_is_reverse(h)) {
                            std::cerr << "[smoothxg::main::create_consensus_graph] error: node " << node_id
                            << " orientation in the consensus graph is " << consensus->get_is_reverse(consensus_h)
                            << " but actually should be " << smoothed.get_is_reverse(h)  << std::endl;
                            exit(1);
                        }
                        // sequence comparison
                        if (consensus->get_sequence(consensus_h) != smoothed.get_sequence(h)) {
                            std::cerr << "[smoothxg::main::create_consensus_graph] error: node " << node_id
                            << " sequence in the consensus graph is " << consensus->get_sequence(consensus_h)
                            << " but actually is " << smoothed.get_sequence(h) << std::endl;
                            exit(1);
                        }
                    });
                }
            });
    */

    consensus->for_each_path_handle(
        [&](const path_handle_t& path) {
            consensus->for_each_step_in_path(path, [&] (const step_handle_t step) {
               if (consensus->has_next_step(step)) {
                   step_handle_t next_step = consensus->get_next_step(step);
                   handle_t h = consensus->get_handle_of_step(step);
                   handle_t next_h = consensus->get_handle_of_step(next_step);
                   if (!consensus->has_edge(h, next_h)) {
                       consensus->create_edge(h, next_h);
                   }
               }
            });
        });

    // unchop the graph
    // this unchop is necessary
    odgi::algorithms::unchop(*consensus, thread_count, false);

    std::cerr << "[smoothxg::create_consensus_graph] removing edges connecting the path with a gap less than " << min_allele_length << "bp" << std::endl;

    // FIXME can this be parallelized?
    // remove edges that are connecting the same path with a gap less than min_allele_length
    // this removes small indel edges (deletions relative to consensus paths)
    ska::flat_hash_set<edge_t> edges_to_remove;
    ska::flat_hash_set<edge_t> edges_to_keep;
    consensus->for_each_path_handle(
        [&](const path_handle_t& path) {
            //std::cerr << "on path " << consensus->get_path_name(path) << std::endl;
            // make a simple positional index
            std::map<std::pair<uint64_t, uint64_t>, int64_t> step_to_pos;
            // check if edges between nodes in the path are jumping more than our scale factor
            int64_t pos = 0;
            consensus->for_each_step_in_path(
                path,
                [&](const step_handle_t& s) {
                    auto h = consensus->get_handle_of_step(s);
                    //std::cerr << "step " << consensus->get_id(h) << ":" << consensus->get_is_reverse(h) << " @ " << pos << " " << as_integers(s)[0] << ":" << as_integers(s)[1] << std::endl;
                    step_to_pos[std::make_pair(as_integers(s)[0],as_integers(s)[1])] = pos;
                    pos += consensus->get_length(consensus->get_handle_of_step(s));
                });
            // now iterate again, but check the edges
            // record which ones jump distances in the path that are < min_allele_length
            //
            consensus->for_each_step_in_path(
                path,
                [&](const step_handle_t& s) {
                    handle_t h = consensus->get_handle_of_step(s);
                    auto key = std::make_pair(as_integers(s)[0],as_integers(s)[1]);
                    int64_t pos = step_to_pos[key] + consensus->get_length(h);
                    // look at all the edges of h
                    // if an edge connects two nodes where one pair of steps is
                    // more distant than min_allele_length
                    // or if the steps are consecutive
                    // we are going to keep it
                    // otherwise, we'll remove it
                    auto in_path
                        = [&](const step_handle_t& step) {
                              return path == consensus->get_path_handle_of_step(step);
                          };
                    if (consensus->get_degree(h, false) > 1) {
                        consensus->follow_edges(
                            h, false,
                            [&](const handle_t& n) {
                                uint64_t c = 0;
                                consensus->for_each_step_on_handle(
                                    n, [&c](const step_handle_t& q) { ++c; });
                                edge_t e = std::make_pair(h, n);
                                if (edges_to_keep.count(e)) {
                                    // ok
                                } else if (c > 1) {
                                    edges_to_keep.insert(e);
                                } else {
                                    // for all steps on the other handle
                                    bool ok = false;
                                    consensus->for_each_step_on_handle(
                                        n,
                                        [&](const step_handle_t& q) {
                                            if (in_path(q)) {
                                                // if the distance is right
                                                auto key = std::make_pair(as_integers(q)[0],
                                                                          as_integers(q)[1]);
                                                int64_t o_pos = step_to_pos[key];
                                                if (o_pos == pos
                                                    || o_pos - pos >= min_allele_length) {
                                                    ok = true;
                                                }
                                            } else {
                                                ok = true;
                                            }
                                        });
                                    edge_t e = std::make_pair(h, n);
                                    if (!ok && !edges_to_keep.count(e)) {
                                        edges_to_remove.insert(e);
                                    } else {
                                        edges_to_remove.erase(e);
                                        edges_to_keep.insert(e);
                                    }
                                }
                            });
                    }
                });
        });

    // remove edges that we selected for removal
    // specifically, these represent indels in one consensus or link
    // shorter than our consensus jump bound
    for (auto& e : edges_to_remove) {
        consensus->destroy_edge(e);
    }

    // force edges in paths
    // TODO this does not seem to be necessary
    /*
    consensus->for_each_path_handle(

        [&](const path_handle_t& path) {
            consensus->for_each_step_in_path(path, [&] (const step_handle_t step) {
               if (consensus->has_next_step(step)) {
                   step_handle_t next_step = consensus->get_next_step(step);
                   handle_t h = consensus->get_handle_of_step(step);
                   handle_t next_h = consensus->get_handle_of_step(next_step);
                   if (!consensus->has_edge(h, next_h)) {
                       consensus->create_edge(h, next_h);
                   }
               }
            });
        });
        */

    odgi::algorithms::unchop(*consensus, thread_count, false);

    auto* copy = new odgi::graph_t();
    graph_deep_copy(consensus, copy);
    delete consensus;
    consensus = copy;

    odgi::algorithms::unchop(*consensus, thread_count, false);

    // remove 0-depth nodes and edges
    auto handles_to_drop = odgi::algorithms::find_handles_exceeding_coverage_limits(*consensus, 1, 0);
    for (auto& handle : handles_to_drop) {
        consensus->destroy_handle(handle);
    }

    odgi::algorithms::unchop(*consensus, thread_count, false);

    uint64_t consensus_nodes = 0;
    uint64_t consensus_length = 0;
    consensus->for_each_handle(
        [&](const handle_t& h) {
            ++consensus_nodes;
            consensus_length += consensus->get_length(h);
        });
    std::cerr << "[smoothxg::create_consensus_graph] final graph length " << consensus_length << "bp " << "in " << consensus_nodes << " nodes" << std::endl;

    return consensus;
}

}
