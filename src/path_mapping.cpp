#include "path_mapping.hpp"

namespace smoothxg {

void path_mappings_t::append(const path_position_range_t& range) {
    mappings->append(range);
}

void path_mappings_t::index(const uint64_t& n_threads) {
    mappings->index(n_threads);
}

uint64_t path_mappings_t::size(void) {
    return mappings->size();
}

path_position_range_t path_mappings_t::get_mapping(const uint64_t& idx) {
    return mappings->read_value(idx);
}

}
