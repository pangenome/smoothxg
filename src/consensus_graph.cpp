#include "consensus_graph.hpp"

namespace smoothxg {

bool operator<(const link_path_t& a,
               const link_path_t& b) {
    auto& a_0 = as_integer(a.from_cons);
    auto& a_1 = as_integer(a.to_cons);
    auto& b_0 = as_integer(b.from_cons);
    auto& b_1 = as_integer(b.to_cons);
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
      << as_integer(a.from_cons) << " "
      << as_integer(a.to_cons) << " "
      << a.length << " "
      << a.hash << " "
      << as_integer(a.path) << " "
      << as_integers(a.begin)[0] << ":" << as_integers(a.begin)[1] << " "
      << as_integers(a.end)[0] << ":" << as_integers(a.end)[1] << ")";
    return o;
}

// prep the graph into a given GFA file
// we'll then build the xg index on top of that in low memory

odgi::graph_t create_consensus_graph(const odgi::graph_t& smoothed,
                                     const std::vector<path_handle_t>& consensus_paths,
                                     const uint64_t& consensus_jump_max,
                                     const uint64_t& thread_count,
                                     const std::string& base) {

    // walk each path
    // record distance since last step on a consensus path
    // record first step handle off a consensus path
    // detect consensus switches, writing the distance to the last consensus step, step
    // into an array of tuples

    
    std::vector<bool> is_consensus(smoothed.get_path_count()+1, false);
    for (auto& path : consensus_paths) {
        is_consensus[as_integer(path)] = true;
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

    std::vector<uint64_t> node_offset;
    node_offset.push_back(0);
    smoothed.for_each_handle(
        [&](const handle_t& h) {
            node_offset.push_back(node_offset.back()+smoothed.get_length(h));
        });

    auto start_in_vector =
        [&](const handle_t& h) {
            if (!smoothed.get_is_reverse(h)) {
                return (int64_t) node_offset[smoothed.get_id(h)-1];
            } else {
                return (int64_t) (node_offset[smoothed.get_id(h)-1]
                                  + smoothed.get_length(h));
            }
        };

    auto end_in_vector =
        [&](const handle_t& h) {
            if (smoothed.get_is_reverse(h)) {
                return (int64_t) node_offset[smoothed.get_id(h)-1];
            } else {
                return (int64_t) (node_offset[smoothed.get_id(h)-1]
                                  + smoothed.get_length(h));
            }
        };

    // consensus path -> consensus path : link_path_t
    mmmulti::set<link_path_t> link_path_ms(base);
    link_path_ms.open_writer();
    
    paryfor::parallel_for<uint64_t>(
        0, non_consensus_paths.size(), thread_count,
        [&](uint64_t idx, int tid) {
            auto& path = non_consensus_paths[idx];
            // for each step in path
            link_path_t link;
            link.path = path;
            path_handle_t last_seen_consensus;
            bool seen_consensus = false;
            smoothed.for_each_step_in_path(
                path,
                [&](const step_handle_t& step) {
                    // check if we're on the step with any consensus
                    handle_t h = smoothed.get_handle_of_step(step);
                    bool on_consensus = false;
                    path_handle_t curr_consensus;
                    smoothed.for_each_step_on_handle(
                        h,
                        [&](const step_handle_t& s) {
                            path_handle_t p = smoothed.get_path_handle_of_step(s);
                            if (is_consensus[as_integer(p)]) {
                                on_consensus = true;
                                curr_consensus = p;
                            }
                        });
                    // if we're on the consensus
                    if (on_consensus) {
                        // we haven't seen any consensus before?
                        if (!seen_consensus) {
                            // we construct the first link path object
                            link.length = 0;
                            link.from_cons = curr_consensus;
                            link.to_cons = curr_consensus;
                            link.begin = step;
                            link.end = step;
                            link.hash = 0;
                            seen_consensus = true;
                            last_seen_consensus = curr_consensus;
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
                            
                            if (link.from_cons == curr_consensus
                                && jump_length < consensus_jump_max) {
                                link.begin = step;
                                link.end = step;
                                link.length = 0;
                            } else { // or it's different
                                // this is when we write a link candidate record
                                link.to_cons = curr_consensus;
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
                                if (as_integer(link.from_cons) > as_integer(link.to_cons)) {
                                    std::swap(link.from_cons, link.to_cons);
                                }
                                link.jump_length = jump_length;
                                link_path_ms.append(link);

                                // reset link
                                link.length = 0;
                                link.from_cons = curr_consensus;
                                link.to_cons = curr_consensus;
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
        });

    link_path_ms.index(thread_count);

    // collect sets of link paths that refer to the same consensus path pairs
    // and pick which one to keep in the consensus graph

    std::vector<link_path_t> consensus_links;
    std::vector<link_path_t> curr_links;

    path_handle_t curr_from_cons;
    path_handle_t curr_to_cons;

    auto novel_sequence_length =
        [&](const link_path_t& link,
            ska::flat_hash_set<uint64_t> seen_nodes) { // by copy
            uint64_t novel_bp = 0;
            for (auto s = link.begin;
                 s != link.end; s = smoothed.get_next_step(s)) {
                handle_t h = smoothed.get_handle_of_step(s);
                uint64_t i = smoothed.get_id(h);
                if (!seen_nodes.count(i)) {
                    novel_bp += smoothed.get_length(h);
                    seen_nodes.insert(i);
                }
            }
            return novel_bp;
        };

    auto mark_seen_nodes =
        [&](const link_path_t& link,
            ska::flat_hash_set<uint64_t>& seen_nodes) { // by ref
            for (auto s = link.begin;
                 s != link.end; s = smoothed.get_next_step(s)) {
                handle_t h = smoothed.get_handle_of_step(s);
                uint64_t i = smoothed.get_id(h);
                if (!seen_nodes.count(i)) {
                    seen_nodes.insert(i);
                }
            }
        };

    auto compute_best_link =
        [&](const std::vector<link_path_t>& links) {
            std::map<uint64_t, uint64_t> hash_counts;
            std::vector<link_path_t> unique_links;
            for (auto& link : links) {
                //std::cerr << link << std::endl;
                auto& c = hash_counts[link.hash];
                if (c == 0) {
                    unique_links.push_back(link);
                }
                ++c;
            }
            std::map<uint64_t, uint64_t> hash_lengths;
            for (auto& link : links) {
                hash_lengths[link.hash] = link.length;
            }
            for (auto& link : links) {
                std::cerr << link << " " << get_path_seq(smoothed.get_next_step(link.begin), link.end) << std::endl;
            }
            for (auto& c : hash_counts) {
                std::cerr << c.first << " -> " << c.second << std::endl;
            }
            uint64_t best_count = 0;
            uint64_t best_hash;
            for (auto& c : hash_counts) {
                if (c.second > best_count) {
                    best_hash = c.first;
                    best_count = c.second;
                }
            }
            std::cerr << "best hash be " << best_hash << std::endl;
            // save the best link path
            link_path_t most_frequent_link;
            for (auto& link : unique_links) {
                if (link.hash == best_hash) {
                    most_frequent_link = link;
                    break;
                }
            }

            // if we have a 0-length link between consensus ends, add it
            handle_t from_end_fwd
                = smoothed.get_handle_of_step(
                    smoothed.path_back(
                        most_frequent_link.from_cons));
            handle_t to_begin_fwd
                = smoothed.get_handle_of_step(
                    smoothed.path_begin(
                        most_frequent_link.to_cons));
            handle_t from_end_rev = smoothed.flip(to_begin_fwd);
            handle_t to_begin_rev = smoothed.flip(from_end_fwd);
            bool has_perfect_link = false;
            link_path_t perfect_link;
            for (auto& link : unique_links) {
                //path_handle_t from_cons;
                //path_handle_t to_cons;
                // are we just stepping from the end of one to the beginning of the other?
                for (step_handle_t s = link.begin;
                     s != link.end; s = smoothed.get_next_step(s)) {
                    step_handle_t next = smoothed.get_next_step(s);
                    if (next == link.end) break;
                    handle_t b = smoothed.get_handle_of_step(s);
                    handle_t e = smoothed.get_handle_of_step(next);
                    if (b == from_end_fwd && e == to_begin_fwd
                        || b == from_end_rev && e == to_begin_rev) {
                        has_perfect_link = true;
                        perfect_link = link;
                        break;
                    }
                    if (has_perfect_link) {
                        break;
                    }
                }
            }

            ska::flat_hash_set<uint64_t> seen_nodes;
            if (has_perfect_link) {
                consensus_links.push_back(perfect_link);
            } else {
                if (most_frequent_link.from_cons != most_frequent_link.to_cons) {
                    consensus_links.push_back(most_frequent_link);
                    mark_seen_nodes(most_frequent_link, seen_nodes);
                }
            }

            // borked

            
            //for (
            // if the consensus sequences are different,
            // what's the shortest hash length
            /*
            if (most_frequent_link.length > 0) {
                std::vector<std::pair<uint64_t, std::pair<uint64_t, uint64_t>>> hash_length_freqs_p;
                for (auto& h : hash_counts) {
                    hash_lengths_p.push_back(std::make_pair(h.second, h.first));
                }
                std::sort(hash_lengths_p.begin(), hash_lengths_p.end());
                //std::cerr << "shortest link is " << hash_lengths_p.front().second << " " << hash_lengths_p.front().first << std::endl;
                auto& shortest_link = hash_lengths_p.front().second;
                if (shortest_link.length < most_frequent_link.length) {
                    consensus_links.push_back(shortest_link);
                    mark_seen_nodes(shortest_link, seen_nodes);
                }
            }
            */
            
            //if (hash_le) {
                //   what is the shortest connection between the ends of the consensus sequences?
            //}

            // keep at least that one
            // and remove any other links shorter than our consensus_jump_max
            // otherwise, for self links, keep all uniques

            // add the shortest
            
            // and most-frequent link (if different)
            // add other links if they add more than consensus_jump_max bp
            // or jump further than that distance
            //std::unordered_set<uint64_t,
            for (auto& link : unique_links) {
                if (link.hash == best_hash) {
                    continue;
                }
                uint64_t novel_bp = novel_sequence_length(link, seen_nodes);
                if (link.jump_length >= consensus_jump_max
                    || novel_bp >= consensus_jump_max) {
                    consensus_links.push_back(link);
                    mark_seen_nodes(link, seen_nodes);
                }
                // compute novel sequence addition
                
            }
        };
    
    // collect edges by node
    // 
    link_path_ms.for_each_value(
        [&](const link_path_t& v) {
            //std::cerr << "on " << v << " with count " << c << std::endl;
            if (curr_links.empty()) {
                curr_from_cons = v.from_cons;
                curr_to_cons = v.to_cons;
                curr_links.push_back(v);
            } else if (curr_from_cons != v.from_cons
                       || curr_to_cons != v.to_cons) {
                // compute the best link in the set of curr_links
                compute_best_link(curr_links);
                // reset the links
                curr_links.clear();
                curr_from_cons = v.from_cons;
                curr_to_cons = v.to_cons;
                curr_links.push_back(v);
            } else {
                curr_links.push_back(v);
            }
        });

    compute_best_link(curr_links);

    link_path_ms.close_reader();
    std::remove(base.c_str());

    // we need to create a copy of the original graph
    // this sounds memory expensive
    //odgi::graph_t consensus_graph; // = smoothed;
    // build an xp index of the smoothed graph
    // iterate through all blocks
    // fetch the consensus path of the given block
    // we can go left or right!
    // for each path hitting the first or last node of the consensus path:
    // 1. get the step
    // 2. get the next step
    // 3. translate to nucleotide position
    // 4. use the path_nuc_range_block_index to find the block id which is corresponding to the range
    // 5. go from block id to consensus path
    // 6. we can get the first and last handle of the consensus path
    // 7. find out the number of nucleotides for hitted path from:
        // a) left cons_path -> left next_cons_path
        // b) left const_path -> right next_cons_path
        // c) right const_path -> left next_cons_path
        // d) right const_path -> right next_const_path
    // 8. choose the shortest path in nucleotides, note start and end step of that path
    // 9. collect such tuples for all paths that we can hit in a consensus path of a current block
    // we might travel to different blocks, so treat them seperately
    // 10. sort by shortest nucleotide and then create a link path from noted start to end step

    // TODO how to deal with loops?
    // TODO we will create each link "twice", what to do?

    // create links according to sorted set
    //std::vector<path_handle_t> all_consensus_paths = consensus_paths;
    //all_consensus_paths.insert(consensus_paths.end(), link_paths.begin(), link_paths.end());

    // create new consensus graph which only has the consensus and link paths in it
    odgi::graph_t consensus;
    // add the consensus paths first
    for (auto& path : consensus_paths) {
        // create the path
        path_handle_t path_cons_graph = consensus.create_path_handle(smoothed.get_path_name(path));
        handle_t cur_handle_in_cons_graph;
        // add the current node first, then add the step
        smoothed.for_each_step_in_path(path,
                                       [&consensus, &smoothed, &path, &cur_handle_in_cons_graph, &path_cons_graph]
                                       (const step_handle_t& step) {
            handle_t h = smoothed.get_handle_of_step(step);
            handle_t next_handle;
            nid_t node_id = smoothed.get_id(h);
            if (!consensus.has_node(node_id)) {
               cur_handle_in_cons_graph = consensus.create_handle(smoothed.get_sequence(h), node_id);
            } else {
               cur_handle_in_cons_graph = consensus.get_handle(node_id);
            }
            bool rev = smoothed.get_is_reverse(h);
            if (rev) {
                consensus.append_step(path_cons_graph, consensus.flip(cur_handle_in_cons_graph));
            } else {
                consensus.append_step(path_cons_graph, cur_handle_in_cons_graph);
            };
        });
    }
    // add link paths and edges not in the consensus paths
    //std::vector<path_handle_t> link_paths;
    for (auto& link : consensus_links) {
        std::cerr << "making " << "Link_" << smoothed.get_path_name(link.from_cons) << "_" << smoothed.get_path_name(link.to_cons) << std::endl;
        // create link paths and paths
        // make the path name
        if (link.length > 0) {
            stringstream s;
            s << "Link_" << smoothed.get_path_name(link.from_cons) << "_" << smoothed.get_path_name(link.to_cons);
            path_handle_t path_cons_graph = consensus.create_path_handle(s.str());
            handle_t cur_handle_in_cons_graph;
            // add the current node first, then add the step
            for (step_handle_t step = smoothed.get_next_step(link.begin);
                 step != link.end;
                 step = smoothed.get_next_step(step)) {
                handle_t h = smoothed.get_handle_of_step(step);
                handle_t next_handle;
                nid_t node_id = smoothed.get_id(h);
                if (!consensus.has_node(node_id)) {
                    cur_handle_in_cons_graph = consensus.create_handle(smoothed.get_sequence(h), node_id);
                } else {
                    cur_handle_in_cons_graph = consensus.get_handle(node_id);
                }
                bool rev = smoothed.get_is_reverse(h);
                if (rev) {
                    consensus.append_step(path_cons_graph, consensus.flip(cur_handle_in_cons_graph));
                } else {
                    consensus.append_step(path_cons_graph, cur_handle_in_cons_graph);
                }
            }
        }
    }

    // finally add the edges
    consensus.for_each_path_handle(
        [&](const path_handle_t& path) {
            consensus.for_each_step_in_path(path, [&] (const step_handle_t step) {
               if (consensus.has_next_step(step)) {
                   step_handle_t next_step = consensus.get_next_step(step);
                   handle_t h = consensus.get_handle_of_step(step);
                   handle_t next_h = consensus.get_handle_of_step(next_step);
                   if (!consensus.has_edge(h, next_h)) {
                       consensus.create_edge(h, next_h);
                   }
               }
            });
        });

    auto link_steps =
        [&](const step_handle_t& a, const step_handle_t& b) {
            handle_t from = smoothed.get_handle_of_step(a);
            handle_t to = smoothed.get_handle_of_step(b);
            consensus.create_edge(consensus.get_handle(smoothed.get_id(from),
                                                       smoothed.get_is_reverse(from)),
                                  consensus.get_handle(smoothed.get_id(to),
                                                       smoothed.get_is_reverse(to)));
        };

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
    // validate consensus graph
    // number of handles
    smoothed.for_each_path_handle(
            [&](const path_handle_t &p) {
                if (is_consensus[as_integer(p)]) {
                    std::string path_name = smoothed.get_path_name(p);
                    if (!consensus.has_path(path_name)) {
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
                        if (!consensus.has_node(node_id)) {
                            std::cerr << "[smoothxg::main::create_consensus_graph] error: node " << node_id
                            << " not present in the consensus graph!" << std::endl;
                            exit(1);
                        }
                        // node orientation comparison
                        handle_t consensus_h = consensus.get_handle(node_id, smoothed.get_is_reverse(h));
                        if (consensus.get_is_reverse(consensus_h) != smoothed.get_is_reverse(h)) {
                            std::cerr << "[smoothxg::main::create_consensus_graph] error: node " << node_id
                            << " orientation in the consensus graph is " << consensus.get_is_reverse(consensus_h)
                            << " but actually should be " << smoothed.get_is_reverse(h)  << std::endl;
                            exit(1);
                        }
                        // sequence comparison
                        if (consensus.get_sequence(consensus_h) != smoothed.get_sequence(h)) {
                            std::cerr << "[smoothxg::main::create_consensus_graph] error: node " << node_id
                            << " sequence in the consensus graph is " << consensus.get_sequence(consensus_h)
                            << " but actually is " << smoothed.get_sequence(h) << std::endl;
                            exit(1);
                        }
                    });
                }
            });

    // unchop the graph
    odgi::algorithms::unchop(consensus);
    return consensus;
}

}