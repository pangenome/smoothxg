#pragma once
#include "mmmultiset.hpp"
#include "odgi/odgi.hpp"
#include "tempfile.hpp"

namespace smoothxg {

using namespace handlegraph;

struct path_position_range_t {
    path_handle_t base_path = as_path_handle(0);   // base path in input graph
    uint64_t start_pos = 0;        // start position of the range
    uint64_t end_pos = 0;          // end position of the range
    step_handle_t start_step = { 0, 0 };  // start step in the base graph
    step_handle_t end_step = { 0, 0 };    // end step in the base graph
    path_handle_t target_path = as_path_handle(0); // target path in smoothed block graph
    uint64_t block_id = 0;  // the block graph id
};

//bool operator<(const path_position_range_t& a, const path_position_range_t& b);
bool operator<(const path_position_range_t& a, const path_position_range_t& b) {
    auto &a_id = as_integer(a.base_path);
    auto &b_id = as_integer(b.base_path);
    return (a_id < b_id || a_id == b_id && a.start_pos < b.start_pos);
}

struct path_mappings_t {
//std::make_unique<mmmulti::set<path_position_range_t>>(_path_mapping_temp);
    std::unique_ptr<mmmulti::set<path_position_range_t>> mappings;
    std::string _mappings_temp;
    void append(const path_position_range_t& pos_range);
    void index(const uint64_t& n_threads);
    uint64_t size(void);
    path_position_range_t get_mapping(const uint64_t& idx);
    path_mappings_t(void) {
        _mappings_temp = temp_file::create();
        mappings = std::make_unique<mmmulti::set<path_position_range_t>>(_mappings_temp);
        mappings->open_writer();
    }
    ~path_mappings_t(void) {
        mappings->close_reader();
        std::remove(_mappings_temp.c_str());
    }
};

}
