#pragma once

#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <vector>
#include "mmmultimap.hpp"
#include "xg.hpp"
#include "flat_hash_map.hpp"
#include "deps/odgi/src/dset64.hpp"

#include "tempfile.hpp"

namespace smoothxg {

using namespace handlegraph;


inline uint64_t path_rank(const step_handle_t& step) {
    return as_integers(step)[0];
}

inline uint64_t step_rank(const step_handle_t& step) {
    return as_integers(step)[1];
}

struct path_range_t {
    step_handle_t begin = { 0, 0 };
    step_handle_t end = { 0, 0 };
    uint64_t length = 0;
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

struct ranked_path_range_t {
    uint64_t rank = 0;

    path_range_t path_range;

    bool operator<(const ranked_path_range_t& pr) const{
        return rank < pr.rank;
    }
    bool operator>(const ranked_path_range_t& pr) const{
        return rank > pr.rank;
    }
    ranked_path_range_t& operator=(const ranked_path_range_t &pr){
        path_range.begin = pr.path_range.begin;
        path_range.end = pr.path_range.end;
        path_range.length = pr.path_range.length;

        rank = pr.rank;

        return *this;
    }
    bool operator!=(const ranked_path_range_t& pr) const{
        return rank != pr.rank;
    }
};

class blockset_t {
private:
    uint64_t _num_blocks = 0;
    mmmulti::map<uint64_t, ranked_path_range_t>* _blocks = nullptr;

    std::string _path_tmp_blocks;

public:
    explicit blockset_t(const std::string& base) {
        _path_tmp_blocks = temp_file::create(base);

        _num_blocks = 0;
        _blocks = new mmmulti::map<uint64_t, ranked_path_range_t>(_path_tmp_blocks, {0});

        _blocks->open_writer();
    }

    ~blockset_t(){
        //_blocks.close_writer();
        //_blocks.close_reader();
        delete _blocks;
        //std::remove(_path_tmp_blocks.c_str()); // The temp_file is deleted automatically
    }

    [[nodiscard]] uint64_t size() const {
        return _num_blocks;
    }

    void add_block(const uint64_t block_id, block_t& block) {
        _num_blocks += 1;

        for (uint64_t rank = 0; rank < block.path_ranges.size(); ++rank) {
            _blocks->append(
                block_id + 1,
                {rank, block.path_ranges[rank]}
            );
        }
    }

    void index(const uint64_t num_threads) {
        _blocks->index(num_threads, _num_blocks);
    }

    [[nodiscard]] block_t get_block(uint64_t block_id) const {
        block_t block;
        for (auto& ranked_path_range : _blocks->values(block_id + 1)){
            block.path_ranges.push_back(ranked_path_range.path_range);
        }
        return block;
    }
};

// find the boundaries of blocks that we can compress with spoa
// assuming a maximum path length within each block
    void smoothable_blocks(
    const xg::XG& graph,
    blockset_t& blockset,
    const uint64_t& max_block_weight,
    const uint64_t& max_block_path_length,
    const uint64_t& max_path_jump,
    const uint64_t& max_edge_jump,
    const bool& order_paths_from_longest,
    int num_threads);

}
