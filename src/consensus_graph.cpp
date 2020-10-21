#include "consensus_graph.hpp"

namespace smoothxg {

bool operator<(const link_path_t& a,
               const link_path_t& b) {
    /*
    auto& a_0 = as_integer(a.from_cons);
    auto& a_1 = as_integer(a.to_cons);
    auto& b_0 = as_integer(b.from_cons);
    auto& b_1 = as_integer(b.to_cons);
    return a_0 < b_0;
    */
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
                    path_handle_t consensus;
                    smoothed.for_each_step_on_handle(
                        h,
                        [&](const step_handle_t& s) {
                            path_handle_t p = smoothed.get_path_handle_of_step(s);
                            if (is_consensus[as_integer(p)]) {
                                on_consensus = true;
                                consensus = p;
                            }
                        });
                    // if we're on the consensus
                    if (on_consensus) {
                        // we haven't seen any consensus before?
                        if (!seen_consensus) {
                            // we construct the first link path object
                            link.length = 0;
                            link.from_cons = consensus;
                            link.to_cons = consensus;
                            link.begin = step;
                            link.end = step;
                            link.hash = 0;
                            seen_consensus = true;
                            last_seen_consensus = consensus;
                        } else {
                            /*
                            if (last_seen_consensus != consensus) {
                                std::cerr << "path " << smoothed.get_path_name(path) << " switched from " << smoothed.get_path_name(last_seen_consensus) << " to " << smoothed.get_path_name(consensus) << std::endl;
                                last_seen_consensus = consensus;
                            }
                            */

                            // we've seen a consensus before, and it's the same
                            if (link.from_cons == consensus) {
                                link.begin = step;
                                link.end = step;
                                link.length = 0;
                            } else { // or it's different
                                // this is when we write a link candidate record
                                link.to_cons = consensus;
                                link.begin = smoothed.get_next_step(link.begin);
                                link.end = step;
                                //std::cerr << "writing to mmset" << std::endl;
                                link.length = get_path_seq_length(link.begin, link.end);
                                link.hash = hash_seq(get_path_seq(link.begin, link.end));
                                link_path_ms.append(link);

                                // reset link
                                link.length = 0;
                                link.from_cons = consensus;
                                link.to_cons = consensus;
                                link.begin = step;
                                link.end = step;
                                link.hash = 0;
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

    std::vector<link_path_t> links;
    std::vector<link_path_t> curr_links;

    path_handle_t curr_from_cons;
    path_handle_t curr_to_cons;

    auto compute_best_link =
        [&](const std::vector<link_path_t>& links) {
            std::map<uint64_t, uint64_t> hash_counts;
            for (auto& link : links) {
                //std::cerr << link << std::endl;
                ++hash_counts[link.hash];
            }
            for (auto& link : links) {
                std::cerr << link << " " << get_path_seq(link.begin, link.end) << std::endl;
            }
            for (auto& c : hash_counts) {
                std::cerr << c.first << " -> " << c.second << std::endl;
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

    std::vector<path_handle_t> link_paths;

    // create links according to sorted set
    std::vector<path_handle_t> all_consensus_paths = consensus_paths;
    all_consensus_paths.insert(consensus_paths.end(), link_paths.begin(), link_paths.end());

    // create new consensus graph which only has the consensus and link paths in it
    odgi::graph_t consensus_graph;
    // add the consensus paths first
    for (auto& path : all_consensus_paths) {
        // create the path
        path_handle_t path_cons_graph = consensus_graph.create_path_handle(smoothed.get_path_name(path));
        handle_t cur_handle_in_cons_graph;
        // add the current node first, then add the step
        smoothed.for_each_step_in_path(path,
                                       [&consensus_graph, &smoothed, &path, &cur_handle_in_cons_graph, &path_cons_graph]
                                       (const step_handle_t& step) {
            handle_t h = smoothed.get_handle_of_step(step);
            handle_t next_handle;
            nid_t node_id = smoothed.get_id(h);
            if (!consensus_graph.has_node(node_id)) {
               cur_handle_in_cons_graph = consensus_graph.create_handle(smoothed.get_sequence(h), node_id);
            } else {
               cur_handle_in_cons_graph = consensus_graph.get_handle(node_id);
            }
            bool rev = smoothed.get_is_reverse(h);
            if (rev) {
                consensus_graph.append_step(path_cons_graph, consensus_graph.flip(cur_handle_in_cons_graph));
            } else {
                consensus_graph.append_step(path_cons_graph, cur_handle_in_cons_graph);
            };
        });
    }
    // finally add the edges
    consensus_graph.for_each_path_handle(
        [&](const path_handle_t& path) {
            consensus_graph.for_each_step_in_path(path, [&] (const step_handle_t step) {
               if (consensus_graph.has_next_step(step)) {
                   step_handle_t next_step = consensus_graph.get_next_step(step);
                   handle_t h = consensus_graph.get_handle_of_step(step);
                   handle_t next_h = consensus_graph.get_handle_of_step(next_step);
                   if (!consensus_graph.has_edge(h, next_h)) {
                       consensus_graph.create_edge(h, next_h);
                   }
               }
            });
        });
    // TODO validate consensus graph
    return consensus_graph;
}

}
