#pragma once

#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <odgi/odgi.hpp>
#include "zstdutil.hpp"

namespace smoothxg {

using namespace handlegraph;

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim);

void graph_deep_copy(odgi::graph_t* source,
                     odgi::graph_t* target);

double handy_parameter(const std::string& value, const double default_value);

std::unique_ptr<odgi::graph_t> get_block_graph(std::vector<std::string *> &block_graphs, const uint64_t &block_id);

void save_block_graph(std::vector<std::string *> &block_graphs, const uint64_t &block_id, const odgi::graph_t *block_graph);

uint64_t modulo(const uint64_t n, const uint64_t d);

}
