#include "consensus_graph.hpp"

namespace smoothxg {

// prep the graph into a given GFA file
// we'll then build the xg index on top of that in low memory

odgi::graph_t create_consensus_graph(ska::flat_hash_map<std::string, IITree<uint64_t , uint64_t>>& happy_tree_friends,
                                     const odgi::graph_t& smoothed,
                                     const std::vector<std::shared_ptr<std::string>>& consensus_names, // pointers to consensus path names for each block
                                     const std::vector<smoothxg::block_t>& blocks,
                                     const std::string& base) {
    // we need to create a copy of the original graph
    // this sounds memory expensive
    odgi::graph_t consensus_graph; // = smoothed;

    // build an xp index of the smoothed graph
    xp::XP path_index;
    path_index.from_handle_graph(smoothed);

    // ...
    
    // iterate through all blocks
    // we can go left or right!
    // for each path, get an end step rank
    // fetch the next step rank and check one of the following:
    // 1. Can we go further into that direction (If we went left, did we reach the leftest step of that path? If we went right, did we reach the rightest step of that path?)
    // 2. Did we land within the same block again? (Is that even possible?)
    // 3. Did we directly hit another block already? -> no link path needed
    // 4. No hit, we continue.
    // As we are progressively traversing all paths step by step "simultaneously", we run into a shortest connection at some point.
    // Then we can finish the loop for this block by adding a link path from the end step of the current path and bock to
    // the start step of the current path and found block
    // TODO Can we enter a block from behind? --> We don't know the exact end step and start step of a block, do we?!
    // TODO SHIT THIS IS NUTS WITHOUT MY ASSUMPTION


    // TODO Do we need to record if two blocks are already connected?

    // TODO Alternatively do a BFS search of the smoothed graph and find out at which steps we hit a block?
    // TODO Then jump ahead the block?
    // TODO Not 100% sure, how one would do this.

    return consensus_graph;
}

}

