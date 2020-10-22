#pragma once

#include <string>
#include <sstream>
#include <unordered_map> // for string hash
#include <odgi/odgi.hpp>
#include <odgi/unchop.hpp>
#include <deps/cgranges/cpp/IITree.h>
#include <deps/odgi/deps/flat_hash_map/flat_hash_map.hpp> // we are using odgi's flat_hash_map here, I suppose this is not a good idea? Maybe rather add the flat_hash_map to smoothxg itself?
#include <mmmultiset.hpp>
#include "paryfor.hpp"
#include "blocks.hpp"

namespace smoothxg {

struct link_path_t {
    path_handle_t from_cons;
    path_handle_t to_cons;
    uint64_t length; // nucleotides
    uint64_t hash;
    step_handle_t begin; // first step off consensus path
    step_handle_t end; // one-past last step
    path_handle_t path;
    bool is_rev; // if we're going forward or reverse on the graph partial order
    uint64_t jump_length; // jump in the partial order
    uint64_t link_rank;
};

bool operator<(const link_path_t& a,
               const link_path_t& b);

ostream& operator<<(ostream& o, const link_path_t& a);


/// build the path_step_rank_ranges -> index_in_blocks_vector
/// flat_hash_map using SKA: KEY: path_name, VALUE: sorted interval_tree using cgranges https://github.com/lh3/cgranges:
/// we collect path_step_rank_ranges and the identifier of an interval is the index of a block in the blocks vector
odgi::graph_t create_consensus_graph(const odgi::graph_t& smoothed,
                                     const std::vector<path_handle_t>& consensus_paths,
                                     const uint64_t& consensus_jump_max,
                                     const uint64_t& thread_count,
                                     const std::string& base);
}
