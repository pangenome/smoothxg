#pragma once

#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <vector>
#include "mmmultimap.hpp"
#include "xg.hpp"

namespace smoothxg {

using namespace handlegraph;

inline uint64_t path_rank(const step_handle_t& step) {
    return as_integers(step)[0];
}

inline uint64_t step_rank(const step_handle_t& step) {
    return as_integers(step)[1];
}

struct path_range_t {
    //uint64_t rank;
    step_handle_t begin;
    step_handle_t end;
    uint64_t length;
    uint64_t nuc_begin;
    uint64_t nuc_end;
};

struct block_t {
    std::vector<handle_t> handles; // hmmm do we need this?
    uint64_t total_path_length = 0; // what of this do we "Really" need?
    //uint64_t max_path_length = 0;
    std::vector<path_range_t> path_ranges;
    bool broken = false;
    bool is_repeat = false;
    //bool is_split = false;
};

class blockset_t {
public:
    mmmulti::map<uint64_t, path_range_t> blocks;
    void add_block(const block_t& block);
    block_t get_block(uint64_t);
    // todo provide comparison operator<() for path_range_t
    //    and, maybe? to sort the mmmultimap properly, comparison for std::pair<uint64_t, path_range_t>
};

// find the boundaries of blocks that we can compress with spoa
// assuming a maximum path length within each block
std::vector<block_t>
smoothable_blocks(
    const xg::XG& graph,
    const uint64_t& max_block_weight,
    const uint64_t& max_path_jump,
    const uint64_t& min_subpath,
    const uint64_t& max_edge_jump,
    const bool& order_paths_from_longest);

}
