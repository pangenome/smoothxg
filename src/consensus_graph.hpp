#pragma once

#include <string>
#include <sstream>
#include <unordered_map> // for string hash
#include <odgi/odgi.hpp>
#include <odgi/unchop.hpp>
#include <deps/cgranges/cpp/IITree.h>
#include <deps/odgi/deps/flat_hash_map/flat_hash_map.hpp> // we are using odgi's flat_hash_map here, I suppose this is not a good idea? Maybe rather add the flat_hash_map to smoothxg itself?
#include <odgi/coverage.hpp>
#include <odgi/groom.hpp>
#include <odgi/gfa_to_handle.hpp>
#include <mmmultiset.hpp>
#include <xg.hpp>
#include "blocks.hpp"
#include "utils.hpp"
#include "atomic_bitvector.hpp"

namespace smoothxg {

struct link_path_t {
    std::string* from_cons_name;
    std::string* to_cons_name;
    path_handle_t from_cons_path;
    path_handle_t to_cons_path;
    uint64_t length; // nucleotides
    uint64_t hash;
    step_handle_t begin; // first step off consensus path
    step_handle_t end; // one-past last step
    path_handle_t path;
    uint64_t jump_length; // jump in the partial order
    uint64_t rank;
};

struct link_range_t {
    nid_t start;
    nid_t end;
    path_handle_t path;
};

struct consensus_spec_t {
    std::string basename;
    int min_allele_len = 0;
    int max_allele_len = 1e6;
    std::string ref_file;
    std::string ref_file_sanitized;
    bool keep_consensus_paths;
    double min_consensus_path_cov = 0;
};

std::string displayname(const consensus_spec_t& spec);

std::vector<consensus_spec_t> parse_consensus_spec(const std::string& spec_str, bool& requires_consensus);

bool operator<(const link_path_t& a,
               const link_path_t& b);

ostream& operator<<(ostream& o, const link_path_t& a);

/// build a consensus graph consisting of consensus paths and link paths between them
odgi::graph_t* create_consensus_graph(const xg::XG &smoothed,
                                      const std::vector<std::string>& consensus_path_names,
                                      const uint64_t& min_allele_length,
                                      const uint64_t& max_allele_length,
                                      const double& min_consensus_path_coverage,
                                      const uint64_t& thread_count,
                                      const std::string& base);
}
