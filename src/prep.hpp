#pragma once

#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <mutex>
#include "ips4o.hpp"
#include "odgi/odgi.hpp"
#include "odgi/chop.hpp"
#include "odgi/topological_sort.hpp"
#include "odgi/gfa_to_handle.hpp"
#include "odgi/dna.hpp"
#include "odgi/xp.hpp"
#include "odgi/path_sgd.hpp"
#include "odgi/threads.hpp"
#include "odgi/groom.hpp"
#include "paryfor.hpp"

namespace smoothxg {

// prep the graph into a given GFA file
// we'll then build the xg index on top of that in low memory

void prep(
    const std::string& gfa_in,
    const std::string& gfa_out,
    const uint64_t& max_node_length = 32);

}
