#pragma once

#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <mutex>
#include "spoa/spoa.hpp"
#include "xg.hpp"
#include "blocks.hpp"
#include "ips4o.hpp"
#include "odgi/odgi.hpp"
#include "odgi/unchop.hpp"
#include "odgi/topological_sort.hpp"
#include "odgi/dna.hpp"
#include "paryfor.hpp"

namespace smoothxg {

struct path_position_range_t {
    path_handle_t base_path; // base path in input graph
    uint64_t start_pos; // start position of the range
    uint64_t end_pos; // end position of the range
    step_handle_t start_step; // start step in the base graph
    step_handle_t end_step; // end step in the base graph
    path_handle_t target_path; // target path in smoothed block graph
    uint64_t target_graph_id; // the block graph id
};

odgi::graph_t smooth(const xg::XG& graph,
                     const block_t& block,
                     std::int8_t poa_m,
                     std::int8_t poa_n,
                     std::int8_t poa_g,
                     std::int8_t poa_e,
                     std::int8_t poa_q,
                     std::int8_t poa_c,
                     const std::string& consensus_name = "");

odgi::graph_t smooth_and_lace(const xg::XG& graph,
                              const std::vector<block_t>& blocks,
                              std::int8_t poa_m,
                              std::int8_t poa_n,
                              std::int8_t poa_g,
                              std::int8_t poa_e,
                              std::int8_t poa_q,
                              std::int8_t poa_c,
                              const std::string& consensus_name = "");

void write_gfa(std::unique_ptr<spoa::Graph>& graph,
               std::ostream& out,
               const std::vector<std::string>& sequence_names,
               bool include_consensus);

void build_odgi(std::unique_ptr<spoa::Graph>& graph,
                odgi::graph_t& output,
                const std::vector<std::string>& sequence_names,
                const std::vector<bool>& aln_is_reverse,
                const std::string& consensus_name,
                bool include_consensus = true);


}
