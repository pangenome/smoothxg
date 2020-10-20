#pragma once

#include <string>
#include <odgi/odgi.hpp>
#include <deps/cgranges/cpp/IITree.h>
#include <deps/odgi/deps/flat_hash_map/flat_hash_map.hpp> // we are using odgi's flat_hash_map here, I suppose this is not a good idea? Maybe rather add the flat_hash_map to smoothxg itself?
#include <odgi/xp.hpp>
#include "blocks.hpp"
#include "odgi/xp.hpp"

namespace smoothxg {

/// build the path_step_rank_ranges -> index_in_blocks_vector
/// flat_hash_map using SKA: KEY: path_name, VALUE: sorted interval_tree using cgranges https://github.com/lh3/cgranges:
/// we collect path_step_rank_ranges and the identifier of an interval is the index of a block in the blocks vector
odgi::graph_t create_consensus_graph(ska::flat_hash_map<std::string, IITree<uint64_t , uint64_t>>& happy_tree_friends,
                                     const odgi::graph_t& smoothed,
                                     const std::vector<std::string>& consensus_names, // pointers to consensus path names for each block
                                     const std::vector<smoothxg::block_t>& blocks,
                                     const std::string& base);
}
