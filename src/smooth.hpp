#pragma once

#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include "spoa/spoa.hpp"
#include "xg.hpp"
#include "blocks.hpp"
#include "odgi/odgi.hpp"
#include "odgi/unchop.hpp"
#include "odgi/topological_sort.hpp"

namespace smoothxg {

odgi::graph_t smooth(const xg::XG& graph,
                     const block_t& block,
                     std::ostream& out,
                     const std::string& consensus_name = "");

void write_gfa(std::unique_ptr<spoa::Graph>& graph,
               std::ostream& out,
               const std::vector<std::string>& sequence_names,
               bool include_consensus);

void build_odgi(std::unique_ptr<spoa::Graph>& graph,
                odgi::graph_t& output,
                const std::vector<std::string>& sequence_names,
                const std::string& consensus_name,
                bool include_consensus = true);


}
