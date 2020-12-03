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
    step_handle_t begin;
    step_handle_t end;
    uint64_t length = 0;
    uint64_t nuc_begin = 0;
    uint64_t nuc_end = 0;

    uint64_t rank = 0;

    bool operator<(const path_range_t& pr) const{
       return rank < pr.rank;
    }
    bool operator>(const path_range_t& pr) const{
        return rank > pr.rank;
    }
    path_range_t& operator=(const path_range_t &pr){
        begin = pr.begin;
        end = pr.end;
        length = pr.length;
        nuc_begin = pr.nuc_begin;
        nuc_end = pr.nuc_end;

        rank = pr.rank;

        return *this;
    }
    bool operator!=(const path_range_t& pr) const{
        return rank != pr.rank;
    }
};

struct block_t {
    //std::vector<handle_t> handles;  // Do we need this? Yes, but not here.
    //uint64_t total_path_length = 0; // Do we need this? Yes, but not here.
    //uint64_t max_path_length = 0; // Not used.
    std::vector<path_range_t> path_ranges;
    //bool broken = false;        // Not used.
    //bool is_repeat = false;     // Not used.
    //bool is_split = false;      // Not used.
};

class blockset_t {
private:
    uint64_t _num_blocks;
    mmmulti::map<uint64_t, path_range_t> _blocks;

    std::string _path_tmp_blocks;
public:
    blockset_t(const std::string& dir_work = "") {
        _path_tmp_blocks = dir_work + "temp.blocks";
        std::remove(_path_tmp_blocks.c_str());

        _num_blocks = 0;
        _blocks.set_base_filename(_path_tmp_blocks);

        _blocks.open_writer();
    }

    ~blockset_t(){
        std::remove(_path_tmp_blocks.c_str());
    }

    [[nodiscard]] uint64_t size() const {
        return _num_blocks;
    }

    void add_block(const uint64_t block_id, block_t& block) {
        _num_blocks += 1;

        for (uint64_t i = 0; i < block.path_ranges.size(); ++i) {
            block.path_ranges[i].rank = i;
            _blocks.append(block_id + 1, block.path_ranges[i]);
        }
    }

    void index(const uint64_t num_threads, const uint64_t max_block_id) {
        _blocks.index(num_threads, max_block_id + 1);
    }

    block_t get_block(uint64_t block_id) {
        block_t block;
        block.path_ranges = _blocks.values(block_id + 1);
        return block;
    }
};

// find the boundaries of blocks that we can compress with spoa
// assuming a maximum path length within each block
    void smoothable_blocks(
    const xg::XG& graph,
    blockset_t& blockset,
    const uint64_t& max_block_weight,
    const uint64_t& max_path_jump,
    const uint64_t& min_subpath,
    const uint64_t& max_edge_jump,
    const bool& order_paths_from_longest);

}
