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
#include "odgi/groom.hpp"
#include "paryfor.hpp"

namespace smoothxg {

void cleanup(
    odgi::graph_t& graph,
    const float& p_sgd_min_term_updates,
    const bool& toposort,
    const uint64_t& num_threads);

}
