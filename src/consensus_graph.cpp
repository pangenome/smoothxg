#include <deps/odgi/src/odgi.hpp>
#include "consensus_graph.hpp"

namespace smoothxg {

bool operator<(const link_path_t& a,
               const link_path_t& b) {
    auto& a_0 = as_integer(a.from_cons_path);
    auto& a_1 = as_integer(a.to_cons_path);
    auto& b_0 = as_integer(b.from_cons_path);
    auto& b_1 = as_integer(b.to_cons_path);
    auto a_0_p = a.from_cons_part;
    auto a_1_p = a.to_cons_part;
    auto b_0_p = b.from_cons_part;
    auto b_1_p = b.to_cons_part;
    return std::tie(a_0, a_0_p, a_1, a_1_p, a.length, a.jump_length, a.begin, a.end, as_integer(a.path))
        < std::tie(b_0, b_0_p, b_1, b_1_p, b.length, b.jump_length, b.begin, b.end, as_integer(b.path));
}

ostream& operator<<(ostream& o, const link_path_t& a) {
    o << "("
      << a.from_cons_name << " "
      << a.to_cons_name << " "
      << as_integer(a.from_cons_path) << " "
      << (int)a.from_cons_part << " "
      << as_integer(a.to_cons_path) << " "
      << (int)a.to_cons_part << " "
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

    auto get_path_seq =
        [&](const step_handle_t& begin,
            const step_handle_t& end) {
            std::string seq;
            for (step_handle_t i = begin; i != end; i = smoothed.get_next_step(i)) {
                seq.append(smoothed.get_sequence(smoothed.get_handle_of_step(i)));
            }
            return seq;
        };

    auto hash_link =
        [&](const link_path_t& a, const std::string& seq) {
            stringstream h;
            h << as_integer(a.from_cons_path) << "_"
              << a.from_cons_part << "_"
              << smoothed.get_id(smoothed.get_handle_of_step(a.begin)) << "_"
              << as_integer(a.to_cons_path) << "_"
              << a.to_cons_part << "_"
              << smoothed.get_id(smoothed.get_handle_of_step(a.end)) << "_"
              << a.length << "_"
              << a.jump_length << "_"
              << seq;
            //std::cerr << h.str() << std::endl;
            //return std::hash(std::tie(as_integer(a.from_cons_path), as_integer(a.to_cons_path), a.begin, a.end, a.length, a.jump_length, seq));
            return std::hash<std::string>{}(h.str());
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

    auto consensus_part =
        [&](const path_handle_t& consensus,
            const handle_t& handle) {
            // find the consensus step on this handle
            int64_t curr_pos = -1;
            smoothed.for_each_step_on_handle(
                handle, [&](const step_handle_t& step) {
                    if (smoothed.get_path_handle_of_step(step) == consensus) {
                        curr_pos = smoothed.get_position_of_step(step);
                    }
                });
            // determine what part of the path we're in
            // the first and last 1/8ths are the "ends"
            // the rest is the middle
            if (curr_pos >= 0) {
                auto cons_len = smoothed.get_path_length(consensus);
                if (curr_pos < cons_len / 8) {
                    return path_part_t::begin;
                } else if (curr_pos < cons_len - (cons_len / 8)) {
                    return path_part_t::middle;
                } else {
                    return path_part_t::end;
                }
            } else {
                return path_part_t::middle;
            }
        };

    auto consensus_distance =
        [&](const path_handle_t& consensus,
            const handle_t& last,
            const handle_t& curr) {
            // get position at end of
            int64_t start_pos = -1;
            int64_t end_pos = -1;
            smoothed.for_each_step_on_handle(
                last, [&](const step_handle_t& step) {
                          if (smoothed.get_path_handle_of_step(step) == consensus) {
                              start_pos = smoothed.get_position_of_step(step);
                              // TODO account for offsets
                              if (last == smoothed.get_handle_of_step(step)) {
                                  start_pos += smoothed.get_length(last);
                              }
                          }
                      });
            smoothed.for_each_step_on_handle(
                curr, [&](const step_handle_t& step) {
                          if (smoothed.get_path_handle_of_step(step) == consensus) {
                              end_pos = smoothed.get_position_of_step(step);
                              if (curr == smoothed.flip(smoothed.get_handle_of_step(step))) {
                                  end_pos += smoothed.get_length(last);
                              }
                          }
                      });
//#pragma omp critical (cerr)
            //std::cerr << "got consensus distance " << std::abs(end_pos - start_pos) << std::endl;
            if (start_pos >= 0 && end_pos >= 0) {
                return std::abs(end_pos - start_pos);
            } else {
                return (int64_t)-1;
            }
        };

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
            // .... currently we are using the "last" one we find (but there's only one)
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
                    link.from_cons_part = consensus_part(curr_consensus, h);
                    link.to_cons_path = curr_consensus;
                    link.to_cons_part = link.from_cons_part;
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
                    auto curr_start_fwd = start_in_vector(curr_handle);
                    auto last_end_fwd = end_in_vector(last_handle);
                    auto curr_start_rev = start_in_vector(smoothed.flip(curr_handle));
                    auto last_end_rev = end_in_vector(smoothed.flip(last_handle));
                    /*
#pragma omp critical (cerr)
                    std::cerr << "consensus_dist " << std::abs(curr_start_fwd - last_end_fwd)
                              << " " << consensus_distance(curr_consensus, last_handle, curr_handle)
                              << std::endl;
                    */
                    uint64_t jump_length = (link.from_cons_path == curr_consensus ?
                                            std::min(std::abs(curr_start_fwd - last_end_fwd),
                                                     // distance between last and current on the given path
                                                     consensus_distance(curr_consensus, last_handle, curr_handle))
                                            : 0);
                                            /*
                                            : std::min(std::min(std::abs(curr_start_fwd - last_end_fwd),
                                                             std::abs(curr_start_rev - last_end_rev)),
                                                    std::min(std::abs(curr_start_fwd - last_end_rev),
                                                             std::abs(curr_start_rev - last_end_fwd))));
                                            */
                    // TODO: don't just look at min_allele_length, but also consider allele frequency
                    if (link.from_cons_path == curr_consensus && jump_length < min_allele_length) {
                        link.begin = step;
                        link.end = step;
                        link.length = 0;
                    } else { // or it's different
                        // this is when we write a link candidate record
                        link.to_cons_name = cons_path_ptr[as_integer(curr_consensus)];
                        link.to_cons_path = curr_consensus;
                        //link.begin = smoothed.get_next_step(link.begin);
                        //link.end = smoothed.get_next_step(step);
                        link.end = step;
                        //std::cerr << "writing to mmset" << std::endl;
                        std::string seq = get_path_seq(smoothed.get_next_step(link.begin), link.end);
                        link.length = seq.length();

                        link.jump_length = jump_length;
                        link.hash = hash_link(link, seq);
                        // TODO flip things around to a canonical orientation
                        // to avoid funny effects from inverting paths
                        //
                        auto h_b = smoothed.get_handle_of_step(link.begin);
                        auto h_e = smoothed.get_handle_of_step(link.end);
                        bool rev_b = smoothed.get_is_reverse(h_b);
                        bool rev_e = smoothed.get_is_reverse(h_e);
                        nid_t id_b = smoothed.get_id(h_b);
                        nid_t id_e = smoothed.get_id(h_e);
                        // part of the consensus path, beginning, middle, end
                        auto part_b = consensus_part(link.from_cons_path, h_b);
                        auto part_e = consensus_part(link.to_cons_path, h_e);
                        if (rev_b && rev_e
                            || ((rev_b || rev_e) && std::tie(id_b,part_b) > std::tie(id_e,part_e))) {
                            std::swap(link.from_cons_path, link.to_cons_path);
                            std::swap(link.from_cons_part, link.to_cons_part);
                            std::swap(link.from_cons_name, link.to_cons_name);
                        }
                        link_path_ms->append(link);
                        is_there_something.store(true);

                        // reset link
                        link.length = 0;
                        link.from_cons_name = cons_path_ptr[as_integer(curr_consensus)];
                        link.to_cons_name = cons_path_ptr[as_integer(curr_consensus)];
                        link.from_cons_path = curr_consensus;
                        link.to_cons_path = curr_consensus;
                        link.from_cons_part = consensus_part(curr_consensus, h);
                        link.to_cons_part = link.from_cons_part;
                        link.begin = step;
                        link.end = step;
                        link.hash = 0;
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
                    novel_bp += !seen_nodes.count(i) ? graph.get_length(h) : 0;
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
                        //seen_nodes.insert(i);
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
                    seen_nodes.insert(i);
                }
            };

        auto mark_seen_node_range =
            [&](const step_handle_t begin,
                const step_handle_t end,
                ska::flat_hash_set<uint64_t> &seen_nodes, // by ref
                const xg::XG &graph) {
                uint64_t min_id = std::numeric_limits<uint64_t>::max();
                uint64_t max_id = std::numeric_limits<uint64_t>::min();
                for (auto s = begin;
                     s != end; s = graph.get_next_step(s)) {
                    handle_t h = graph.get_handle_of_step(s);
                    uint64_t i = graph.get_id(h);
                    min_id = std::min(i, min_id);
                    max_id = std::max(i, max_id);
                }
                for (uint64_t i = min_id; i <= max_id; ++i) {
                    seen_nodes.insert(i);
                }
            };
        // collect sets of link paths that refer to the same consensus path pairs
        // and pick which one to keep in the consensus graph

        std::vector<std::vector<link_path_t>> thread_consensus_links(thread_count);
        std::vector<std::vector<std::pair<handle_t, handle_t>>> thread_perfect_edges(thread_count);
        auto compute_link_paths =
                [&](const std::vector<link_path_t> &links) {
                    uint64_t tid = omp_get_thread_num();
                    std::map<uint64_t, uint64_t> hash_counts;
                    std::vector<link_path_t> unique_links;
                    for (auto &link : links) {
                        auto &c = hash_counts[link.hash];
                        if (c == 0) {
                            unique_links.push_back(link);
                        }
                        ++c;
                    }
                    std::sort(unique_links.begin(),
                              unique_links.end(),
                              [&](const link_path_t& a,
                                  const link_path_t& b) {
                                  auto& a_0 = as_integer(a.from_cons_path);
                                  auto& a_1 = as_integer(a.to_cons_path);
                                  auto& b_0 = as_integer(b.from_cons_path);
                                  auto& b_1 = as_integer(b.to_cons_path);
                                  auto a_0_p = a.from_cons_part;
                                  auto a_1_p = a.to_cons_part;
                                  auto b_0_p = b.from_cons_part;
                                  auto b_1_p = b.to_cons_part;
                                  return std::tie(hash_counts[a.hash], a_0, a_0_p, a_1, a_1_p, a.begin, a.end, a.length, a.jump_length, as_integer(a.path))
                                      > std::tie(hash_counts[b.hash], b_0, b_0_p, b_1, b_1_p, b.begin, b.end, b.length, b.jump_length, as_integer(b.path));
                              });
                    // save the best link path
                    link_path_t most_frequent_link = unique_links.front();

                    // keep all edges that directly hop into another consensus
                    uint64_t perfect_edge_count = 0;
                    ska::flat_hash_set<uint64_t> seen_nodes;

                    auto link_cons_end =
                        [&](const step_handle_t& s_end, bool go_rev, const path_handle_t& target_path) {
                            auto cons_end = smoothed.get_handle_of_step(s_end);
                            smoothed.follow_edges(
                                cons_end, go_rev,
                                [&](const handle_t& n) {
                                    smoothed.for_each_step_on_handle(
                                        n,
                                        [&](const step_handle_t& s) {
                                            if (smoothed.get_path_handle_of_step(s)
                                                == target_path) {
                                                auto p = (!go_rev ? std::make_pair(cons_end, n) :
                                                          std::make_pair(n, cons_end));
                                                thread_perfect_edges[tid].push_back(p);
                                                //mark_seen_nodes(, most_frequent_link.end, seen_nodes, smoothed);
                                                uint64_t i = smoothed.get_id(cons_end);
                                                if (!seen_nodes.count(i)) { seen_nodes.insert(i); }
                                                i = smoothed.get_id(n);
                                                if (!seen_nodes.count(i)) { seen_nodes.insert(i); }
                                                ++perfect_edge_count;
                                            }
                                        });
                                });
                        };

                    auto a = std::make_pair(most_frequent_link.from_cons_path, most_frequent_link.from_cons_part);
                    auto b = std::make_pair(most_frequent_link.to_cons_path, most_frequent_link.to_cons_part);
                    bool diff_consensi = a != b;
                    //std::cerr << "getting links for " << smoothed.get_path_name(a) << " and " << smoothed.get_path_name(b) << std::endl;

                    if (diff_consensi) {
                        //handle_t from_end_fwd
                        link_cons_end(smoothed.path_back(a.first), false, b.first);
                        link_cons_end(smoothed.path_begin(a.first), true, b.first);
                        link_cons_end(smoothed.path_back(b.first), false, a.first);
                        link_cons_end(smoothed.path_begin(b.first), true, a.first);
                        
                        link_cons_end(smoothed.path_back(a.first), true, b.first);
                        link_cons_end(smoothed.path_begin(a.first), false, b.first);
                        link_cons_end(smoothed.path_back(b.first), true, a.first);
                        link_cons_end(smoothed.path_begin(b.first), false, a.first);
                    }

                    mark_seen_nodes(smoothed.path_begin(a.first), smoothed.path_end(a.first), seen_nodes, smoothed);
                    mark_seen_nodes(smoothed.path_begin(b.first), smoothed.path_end(b.first), seen_nodes, smoothed);

                    auto& save_links = thread_consensus_links[tid];
                    uint64_t link_rank = 0;

                    // this part attempts to preserve connectivity between consensus sequences
                    // we're preserving the consensus graph topology
                    {
                        // todo iterate through each pair of start/end positions
                        // and keep only the best
                        //std::cerr << "got " << perfect_edge_count << " perfect edges and " << unique_links.size() << " unique links" << std::endl;
                        uint64_t links_to_add = std::min((uint64_t)std::max((int64_t)0, (int64_t)2-(int64_t)perfect_edge_count),
                                                         (uint64_t)unique_links.size());
                        uint64_t zero_links = 0;
                        uint64_t regular_links = 0;
                        for (uint64_t i = 0; i < links_to_add; ++i) {
                            auto& link = unique_links[i];
                            if (link.length == 0) {
                                if (diff_consensi) {
                                    thread_perfect_edges[tid].push_back(
                                        std::make_pair(
                                            smoothed.get_handle_of_step(link.begin),
                                            smoothed.get_handle_of_step(link.end)));
                                    ++zero_links;
                                }
                            } else {
                                link.rank = link_rank++;
                                save_links.push_back(link);
                                mark_seen_nodes(link.begin, link.end, seen_nodes, smoothed);
                                ++regular_links;
                            }
                        }
                        // this part collects sequences that diverge from a consensus for the
                        // min_allele_length which is the "variant scale factor" of the algorithm
                        // this preserves novel non-consensus sequences greater than this length
                        // TODO: this could be made to respect allele frequency
                        for (uint64_t i = links_to_add; i < unique_links.size(); ++i) {
                            auto& link = unique_links[i];
                            uint64_t largest_novel_gap_bp = largest_novel_gap(link.begin, link.end, seen_nodes, smoothed);
                            uint64_t novel_bp = novel_sequence_length(link.begin, link.end, seen_nodes, smoothed);
                            //uint64_t step_count = get_step_count(link.begin, link.end, smoothed);
                            // this complex filter attempts to keep representative link paths for indels above our min_allele_length
                            // we either need the jump length (measured in terms of delta in our graph vector) to be over our jump max
                            // *and* the link path should be empty or mostly novel
                            // *or* we're adding in the specified amount of novel_bp of sequence
                            if (link.length == 0 && diff_consensi) {
                                thread_perfect_edges[tid].push_back(
                                    std::make_pair(
                                        smoothed.get_handle_of_step(link.begin),
                                        smoothed.get_handle_of_step(link.end)));
                                ++zero_links;
                            } else if (link.length == novel_bp
                                && ((most_frequent_link.from_cons_path == most_frequent_link.to_cons_path
                                     && link.jump_length >= min_allele_length
                                     && link.jump_length < max_allele_length
                                     && (link.length == 0
                                         || (novel_bp >= min_allele_length
                                             && largest_novel_gap_bp == novel_bp)))
                                    || (novel_bp >= min_allele_length
                                        && novel_bp < max_allele_length))) {
                                link.rank = link_rank++;
                                save_links.push_back(link);
                                mark_seen_nodes(link.begin, link.end, seen_nodes, smoothed);
                                ++regular_links;
                            }
                        }
                        //std::cerr << "got 0-length links " << zero_links << " and " << regular_links << " regulars" << std::endl;
                    }
                };

        std::cerr << "[smoothxg::create_consensus_graph] finding consensus links" << std::endl;
        // collect edges by noder
        // into groups that we will evaluate in parallel
        std::vector<std::pair<uint64_t, uint64_t>> link_groups;
        std::pair<link_path_t,uint64_t> last = std::make_pair(link_path_ms->read_value(0), 0);
        for (size_t i = 1; i < link_path_ms->size(); ++i) {
            link_path_t curr = link_path_ms->read_value(i);
            if (last.first.from_cons_path != curr.from_cons_path
                || last.first.from_cons_part != curr.from_cons_part
                || last.first.to_cons_path != curr.to_cons_path
                || last.first.to_cons_part != curr.to_cons_part) {
                link_groups.push_back(std::make_pair(last.second, i));
                last = std::make_pair(curr, i);
            }
        }
        link_groups.push_back(std::make_pair(last.second, link_path_ms->size()));
        // run the parallel computation of link paths
#pragma omp parallel for schedule(dynamic,1)
        for (auto& group : link_groups) {
            std::vector<link_path_t> curr_links;
            auto i = group.first;
            auto j = group.second;
            for ( ; i != j ; ++i) {
                curr_links.push_back(link_path_ms->read_value(i));
            }
            compute_link_paths(curr_links);
        }
        for (auto& cons_links : thread_consensus_links) {
            consensus_links.reserve(consensus_links.size() + distance(cons_links.begin(), cons_links.end()));
            consensus_links.insert(consensus_links.end(), cons_links.begin(), cons_links.end());
        }
        for (auto& perf_edges : thread_perfect_edges) {
            perfect_edges.reserve(perfect_edges.size() + distance(perf_edges.begin(), perf_edges.end()));
            perfect_edges.insert(perfect_edges.end(), perf_edges.begin(), perf_edges.end());
        }
    }

    link_path_ms->close_reader();
    link_path_ms.reset();
    std::remove(base_mmset.c_str());

    ips4o::parallel::sort(consensus_links.begin(), consensus_links.end());

    // create new consensus graph which only has the consensus and link paths in it
    auto* consensus = new odgi::graph_t();
    consensus->set_number_of_threads(thread_count);

    std::vector<link_path_t> keep_consensus_links;
    std::vector<std::string> link_path_names;

    // add the consensus paths first
    // FIXME could this be run in parallel? --> create the path in an extra for
    {
        std::vector<bool> link_seen_nodes(smoothed.get_node_count()+1, false);

        std::cerr << "[smoothxg::create_consensus_graph] adding consensus paths" << std::endl;
        for (auto& path : consensus_paths) {
            // create the path
            path_handle_t path_cons_graph = consensus->create_path_handle(smoothed.get_path_name(path));
            handle_t cur_handle_in_cons_graph;
            // add the current node first, then add the step
            smoothed.for_each_step_in_path(path,
                                           [&consensus, &smoothed, &cur_handle_in_cons_graph, &path_cons_graph, &link_seen_nodes]
                                           (const step_handle_t& step) {
                                               handle_t h = smoothed.get_handle_of_step(step);
                                               bool rev = smoothed.get_is_reverse(h);
                                               nid_t node_id = smoothed.get_id(h);
                                               link_seen_nodes[node_id] = true;
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

        auto mark_seen_nodes =
            [&](const step_handle_t begin,
                const step_handle_t end,
                std::vector<bool> &seen_nodes, // by ref
                const xg::XG &graph) {
                for (auto s = begin;
                     s != end; s = graph.get_next_step(s)) {
                    handle_t h = graph.get_handle_of_step(s);
                    uint64_t i = graph.get_id(h);
                    seen_nodes[i] = true;
                }
            };

        auto largest_novel_gap =
            [&](const step_handle_t begin,
                const step_handle_t end,
                const std::vector<bool>& seen_nodes,
                const xg::XG &graph) {
                uint64_t novel_bp = 0;
                uint64_t largest_gap = 0;
                for (auto s = begin;
                     s != end; s = graph.get_next_step(s)) {
                    handle_t h = graph.get_handle_of_step(s);
                    uint64_t i = graph.get_id(h);
                    if (!seen_nodes[i]) {
                        novel_bp += graph.get_length(h);
                    } else {
                        largest_gap = std::max(novel_bp, largest_gap);
                        novel_bp = 0;
                    }
                }
                return largest_gap;
            };

        auto novel_sequence_length =
            [&](const step_handle_t begin,
                const step_handle_t end,
                std::vector<bool> &seen_nodes, // by ref
                const xg::XG &graph) {
                uint64_t novel_bp = 0;
                for (auto s = begin;
                     s != end; s = graph.get_next_step(s)) {
                    handle_t h = graph.get_handle_of_step(s);
                    uint64_t i = graph.get_id(h);
                    novel_bp += !seen_nodes[i] ? graph.get_length(h) : 0;
                }
                return novel_bp;
            };

        auto novel_start =
            [&](const step_handle_t begin,
                const step_handle_t end,
                std::vector<bool> &seen_nodes, // by ref
                const xg::XG &graph) {
                auto s = begin;
                for ( ; s != end; s = graph.get_next_step(s)) {
                    handle_t h = graph.get_handle_of_step(s);
                    uint64_t i = graph.get_id(h);
                    if (seen_nodes[i]) break;
                }
                return std::make_pair(begin, s);
            };

        auto novel_end =
            [&](const step_handle_t begin,
                const step_handle_t end,
                std::vector<bool> &seen_nodes, // by ref
                const xg::XG &graph) {
                auto s = graph.get_previous_step(end);
                for ( ; ; s = graph.get_previous_step(s)) {
                    handle_t h = graph.get_handle_of_step(s);
                    uint64_t i = graph.get_id(h);
                    if (seen_nodes[i]) break;
                    if (s == begin) break;
                }
                return std::make_pair(s, end);
            };

        std::cerr << "[smoothxg::create_consensus_graph] adding link paths: adding paths and nodes" << std::endl;

        // add link paths and edges not in the consensus paths
        for (auto& link : consensus_links) {
            //std::cerr << "making " << "Link_" << smoothed.get_path_name(link.from_cons_path) << "_" << smoothed.get_path_name(link.to_cons_path) << "_" << link.rank << std::endl;
            // create link paths and paths
            if (link.length > 0) {
                auto& novel_link = link;
                // check if this link path is completely novel, to avoid adding "accordions" and SNP arrays to the consensus
                uint64_t novel_bp = novel_sequence_length(link.begin, link.end, link_seen_nodes, smoothed);
                //std::cerr << "link length vs novel bp " << link.length << " " << novel_bp << std::endl;
                std::vector<link_path_t> links;
                if (link.length > novel_bp) {
                    uint64_t largest_novel_gap_bp = largest_novel_gap(link.begin, link.end, link_seen_nodes, smoothed);
                    // attemp to add the link if it's very novel, but only after we trim it into pieces
                    if (largest_novel_gap_bp > 0.95 * link.length) {
                        // here we need to convert the link path into two paths that continue until they're not novel anymore
                        //std::cerr << "not novel enough" << std::endl;
                        auto start_novel = novel_start(link.begin, link.end, link_seen_nodes, smoothed);
                        auto end_novel = novel_end(link.begin, link.end, link_seen_nodes, smoothed);
                        if (start_novel.first != start_novel.second) {
                            auto link_start = link;
                            link_start.begin = start_novel.first;
                            link_start.end = start_novel.second;
                            links.push_back(link_start);
                            //std::cerr << "adding start link" << std::endl;
                        }
                        if (end_novel.first != end_novel.second) {
                            auto link_end = link;
                            link_end.begin = end_novel.first;
                            link_end.end = end_novel.second;
                            links.push_back(link_end);
                            //std::cerr << "adding end link" << std::endl;
                        }
                    }
                } else {
                    links.push_back(link);
                }
                uint64_t i = 0;
                for (auto& link : links) {
                    mark_seen_nodes(link.begin, link.end, link_seen_nodes, smoothed);
                    // make the path name
                    stringstream s;
                    s << "Link_" << *novel_link.from_cons_name << "_"
                      << (char)novel_link.from_cons_part << "_"
                      << *novel_link.to_cons_name << "_"
                      << (char)novel_link.to_cons_part << "_"
                      << novel_link.rank << "_" << i++;
                    //std::cerr << "adding link " << s.str() << std::endl;
                    assert(!consensus->has_path(s.str()));
                    path_handle_t path_cons_graph = consensus->create_path_handle(s.str());
                    //link_paths.push_back(path_cons_graph);

                    // add the current node first, then add the step
                    handle_t cur_handle_in_cons_graph;
                    uint64_t step_count = 0;
                    for (step_handle_t step = (smoothed.has_next_step(novel_link.begin)
                                               ? smoothed.get_next_step(novel_link.begin)
                                               : novel_link.begin);
                         step != novel_link.end;
                         step = smoothed.get_next_step(step)) {
                        handle_t h = smoothed.get_handle_of_step(step);
                        nid_t node_id = smoothed.get_id(h);
                        if (!consensus->has_node(node_id)) {
                            //assert(false);
                            // this must be kept
                            cur_handle_in_cons_graph = consensus->create_handle(smoothed.get_sequence(smoothed.get_handle(node_id)), node_id);
                        } else {
                            cur_handle_in_cons_graph = consensus->get_handle(node_id);
                        }

                        if (smoothed.get_is_reverse(h)) {
                            consensus->append_step(path_cons_graph, consensus->flip(cur_handle_in_cons_graph));
                        } else {
                            consensus->append_step(path_cons_graph, cur_handle_in_cons_graph);
                        }
                        ++step_count;
                    }
                    if (step_count) {
                        link_path_names.push_back(s.str());
                    }
                }
            }
        }
    }

    std::cerr << "[smoothxg::create_consensus_graph] adding link paths: adding edges" << std::endl;
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

    std::cerr << "[smoothxg::create_consensus_graph] adding link paths: adding " << perfect_edges.size() << " perfect edges" << std::endl;
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


    std::cerr << "[smoothxg::create_consensus_graph] adding link paths: walking the paths" << std::endl;
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

    // This is necessary
    odgi::algorithms::unchop(*consensus, thread_count, true);

    auto* copy = new odgi::graph_t();
    graph_deep_copy(consensus, copy);
    delete consensus;
    consensus = copy;

    // remove 0-depth nodes and edges
    std::vector<handle_t> handles_to_drop = odgi::algorithms::find_handles_exceeding_depth_limits(*consensus, 1, 0);
    if (!handles_to_drop.empty()) {
        for (auto& handle : handles_to_drop) {
            consensus->destroy_handle(handle);
        }

        odgi::algorithms::unchop(*consensus, thread_count, true);
    }

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
