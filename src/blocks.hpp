#pragma once

#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <vector>
//#include <IITree.h> // cgranges
//#include "atomic_bitvector.hpp"
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
    step_handle_t begin;
    step_handle_t end;
    uint64_t length;
};

struct block_t {
    std::vector<handle_t> handles;
    uint64_t total_path_length = 0;
    uint64_t max_path_length = 0;
    std::vector<path_range_t> path_ranges;
    bool broken = false;
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
