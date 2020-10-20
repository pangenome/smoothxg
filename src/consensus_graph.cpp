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

    // consensus path -> consensus path : link_path_t
    mmmulti::set<link_path_t> link_path_ms(base);
    link_path_ms.open_writer();
    
    paryfor::parallel_for<uint64_t>(
        0, non_consensus_paths.size(), thread_count,
        [&](uint64_t idx, int tid) {
            auto& path = non_consensus_paths[idx];
            // for each step in path
            link_path_t link;
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
                        } else {
                            // we've seen a consensus before, and it's the same
                            if (link.to_cons == consensus) {
                                // TODO is this the correct place to increment
                                link.length += smoothed.get_length(h);
                            } else { // or it's different
                                // this is when we write a link candidate record
                                link.to_cons = consensus;
                                link.end = step;
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
                    }
                });
        });

    link_path_ms.index(thread_count);

    // create links according to sorted set
    std::vector<path_handle_t> link_paths;
    // TODO

    // create new consensus graph which only has the consensus and link paths in it
    odgi::graph_t consensus_graph;
    // add the consensus paths first
    for (auto& path : consensus_paths) {
        // create the path
        consensus_graph.create_path_handle(smoothed.get_path_name(path));
        // add the current node first, then add the step
        smoothed.for_each_step_in_path(path, [&consensus_graph, &smoothed, &path] (const step_handle_t& step) {
            handle_t h = smoothed.get_handle_of_step(step);
            nid_t node_id = smoothed.get_id(h);
            if (!consensus_graph.has_node(node_id)) {
               consensus_graph.create_handle(smoothed.get_sequence(h), node_id);
            }// TODO somehow insert or set the step ?!
            // consensus_graph.set_step(step, h);
        });
    }
    link_path_ms.close_writer();
    return consensus_graph;
}

}
