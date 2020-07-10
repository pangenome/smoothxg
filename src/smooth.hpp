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
//#include "odgi/unchop.hpp"

namespace smoothxg {

//using namespace handlegraph;
//using nid_t = handlegraph::nid_t;

void smooth(const xg::XG& graph,
            const block_t& block,
            std::ostream& out);

void write_gfa(std::unique_ptr<spoa::Graph>& graph,
               std::ostream& out,
               const std::vector<std::string>& sequence_names,
               bool include_consensus);

void build_odgi(std::unique_ptr<spoa::Graph>& graph,
                odgi::graph_t& output,
                const std::vector<std::string>& sequence_names,
                bool include_consensus);

}
