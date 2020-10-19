#pragma once

#include <string>
#include <deps/cgranges/cpp/IITree.h>
#include <deps/odgi/deps/flat_hash_map/flat_hash_map.hpp> // we are using odgi's flat_hash_map here, I suppose this is not a good idea? Maybe rather add the flat_hash_map to smoothxg itself?
#include "blocks.hpp"
#include "xg.hpp"

namespace smoothxg {

    /// build the path_step_rank_ranges -> path_range_t*
    /// flat_hash_map using SKA: KEY: path_name, VALUE: sorted interval_tree using cgranges https://github.com/lh3/cgranges:
    /// we collect path_step_rank_ranges and the identifier of an interval is a pointer to the corresponding path range type path_range_t*
ska::flat_hash_map<std::string, IITree<uint64_t , uint64_t>> generate_path_nuc_range_block_index(
    std::vector<smoothxg::block_t>& blocks, xg::XG& graph
    );
}