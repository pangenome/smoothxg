#pragma once

#include "blocks.hpp"
#include "deps/abPOA/include/abpoa.h"
#include "deps/abPOA/src/kdq.h"
#include "deps/abPOA/src/utils.h"
#include "ips4o.hpp"
#include "odgi/dna.hpp"
#include "odgi/odgi.hpp"
#include "odgi/topological_sort.hpp"
#include "odgi/unchop.hpp"
#include "paryfor.hpp"
#include "spoa/spoa.hpp"
#include "xg.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <mutex>
#include <sstream>
#include <vector>

namespace smoothxg {

struct path_position_range_t {
    path_handle_t base_path;   // base path in input graph
    uint64_t start_pos;        // start position of the range
    uint64_t end_pos;          // end position of the range
    step_handle_t start_step;  // start step in the base graph
    step_handle_t end_step;    // end step in the base graph
    path_handle_t target_path; // target path in smoothed block graph
    uint64_t target_graph_id;  // the block graph id
};

void write_tsv_for_block(const xg::XG &graph,
                         const block_t &block,
                         const uint64_t &block_id,
                         const std::vector<std::string>& seqs,
                         const std::vector<std::string>& names);

odgi::graph_t smooth_abpoa(const xg::XG &graph, const block_t &block, const uint64_t &block_id,
                           int poa_m, int poa_n, int poa_g,
                           int poa_e, int poa_q, int poa_c,
                           bool local_alignment,
                           std::string *maf,
                           bool banded_alignment,
                           const std::string &consensus_name = "",
                           bool save_block_fastas = false);

odgi::graph_t smooth_spoa(const xg::XG &graph, const block_t &block,
                          const uint64_t &block_id,
                          std::int8_t poa_m, std::int8_t poa_n, std::int8_t poa_g,
                          std::int8_t poa_e, std::int8_t poa_q, std::int8_t poa_c,
                          bool local_alignment,
                          std::string *maf,
                          const std::string &consensus_name = "",
                          bool save_block_fastas = false);

odgi::graph_t smooth_and_lace(const xg::XG &graph,
                              const std::vector<block_t> &blocks,
                              int poa_m, int poa_n,
                              int poa_g, int poa_e,
                              int poa_q, int poa_c,
                              bool local_alignment,
                              int n_threads,
                              std::string &path_output_maf, std::string &maf_header,
                              bool merge_blocks, bool preserve_unmerged_consensus, double contiguous_path_jaccard,
                              bool use_abpoa,
                              const std::string &consensus_name,
                              std::vector<std::string>& consensus_path_names,
                              bool write_fasta_blocks);

void write_gfa(std::unique_ptr<spoa::Graph> &graph, std::ostream &out,
               const std::vector<std::string> &sequence_names,
               bool include_consensus);

void build_odgi(std::unique_ptr<spoa::Graph> &graph, odgi::graph_t &output,
                const std::vector<std::string> &sequence_names,
                const std::vector<bool> &aln_is_reverse,
                const std::string &consensus_name,
                bool include_consensus = true);

void build_odgi_abPOA(abpoa_t *ab, abpoa_para_t *abpt, odgi::graph_t &output,
                      const std::vector<std::string> &sequence_names,
                      const std::vector<bool> &aln_is_reverse,
                      const std::string &consensus_name,
                      bool include_consensus = true);
} // namespace smoothxg
