#include "consensus_graph.hpp"

namespace smoothxg {

// prep the graph into a given GFA file
// we'll then build the xg index on top of that in low memory

odgi::graph_t create_consensus_graph(const odgi::graph_t& smoothed,
                                     const std::vector<path_handle_t>& consensus_paths,
                                     const std::string& base) {

    // walk each path
    // record distance since last step on a consensus path
    // record first step handle off a consensus path
    // detect consensus switches, writing the distance to the last consensus step, step
    // into an array of tuples

    for (auto& path : consensus_paths) {
        std::cerr << "using consensus path " << smoothed.get_path_name(path) << std::endl;
    }
    
    // we need to create a copy of the original graph
    // this sounds memory expensive
    odgi::graph_t consensus_graph; // = smoothed;
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

    // 11. create new consensus graph which only has the consensus and link paths in it

    return consensus_graph;
}

}

