#pragma once

#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <odgi/odgi.hpp>

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

    void graph_deep_copy(odgi::graph_t *source,
                         odgi::graph_t *target);

    std::string* generate_random_string();
}
