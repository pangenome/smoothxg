/** \file smoothxg
 *
 * smooth a graph
 */


#include <omp.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <deps/odgi/src/odgi.hpp>
#include "args.hxx"
#include "sdsl/bit_vectors.hpp"
#include "chain.hpp"
#include "blocks.hpp"
#include "smooth.hpp"
#include "xg.hpp"
#include "prep.hpp"
#include "cleanup.hpp"
#include "breaks.hpp"
#include "utils.hpp"
#include "odgi/odgi.hpp"
#include "consensus_graph.hpp"
#include "rkmh.hpp"
#include <chrono>

using namespace std;
using namespace xg;


int main(int argc, char **argv) {
    args::ArgumentParser parser("smoothxg: collinear block finder and graph consensus generator");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> gfa_in(parser, "FILE", "index the graph in this GFA file", {'g', "gfa-in"});
    args::ValueFlag<std::string> xg_in(parser, "FILE", "read the xg index from this file", {'i', "in"});
    args::ValueFlag<std::string> smoothed_out(parser, "FILE",
                                              "write GFA to this file (not /dev/stdout if consensus graph is made)",
                                              {'o', "smoothed-out"});
    args::ValueFlag<std::string> _smoothed_in_gfa(parser, "FILE", "read GFA from this file", {'F', "smoothed-in"});
    args::ValueFlag<std::string> write_msa_in_maf_format(parser, "FILE",
                                                         "write the multiple sequence alignments (MSAs) in MAF format in this file",
                                                         {'m', "write-msa-in-maf-format"});
    args::ValueFlag<std::string> _ref_paths(parser, "FILE",
                                            "a file listing (one per line) sequences to preserved as paths in the consensus output graphs",
                                            {'P', "ref-paths"});
    //args::Flag _only_ref_paths(parser, "", "use only the reference paths in the consensus graph, ignoring any other consensus paths", {'O', "only-ref-paths"});
    args::ValueFlag<std::string> _write_consensus_path_names(parser, "FILE",
                                                             "write the consensus path names to this file",
                                                             {'f', "write-consensus-path-names"});
    args::ValueFlag<std::string> _read_consensus_path_names(parser, "FILE",
                                                            "don't smooth, just generate the consensus, taking the consensus path names from this file",
                                                            {'H', "consensus-from"});
    //args::ValueFlag<std::string> write_consensus_graph(parser, "BASENAME", "write the consensus graph to BASENAME.cons_[spec].gfa", {'s', "write-consensus-graph"});
    args::ValueFlag<std::string> _consensus_spec(parser, "BASENAME[,min_len[:refs[:(y|n)[:min_cov[:max_len]?]?]?]?]*",
                                                 "consensus graph specification: write the consensus graph to BASENAME.cons_[spec].gfa; where each spec contains at least a min_len parameter (which defines the length of divergences from consensus paths to preserve in the output), optionally a file containing reference paths to preserve in the output, a flag (y/n) indicating whether we should also use the POA consensus paths, a minimum coverage of consensus paths to retain (min_cov), and a maximum allele length (max_len, defaults to 1e6); example: cons,100,1000:refs1.txt:n,1000:refs2.txt:y:2.3:1000000,10000 [default: unset]",
                                                 {'C', "consensus-spec"});

    args::ValueFlag<std::string> _consensus_path_prefix(parser, "PREFIX",
                                                        "prepend the consensus path names with PREFIX [default: Consensus]",
                                                        {'Q', "consensus-prefix"});
    // Merge blocks (for merging MAF blocks and consensus sequences)
    args::Flag merge_blocks(parser, "bool",
                            "merge contiguous MAF blocks in the MAF output and consensus sequences in the smoothed graph",
                            {'M', "merge-blocks"});
    args::Flag vanish_consensus(parser, "bool",
                                "remove the consensus paths from the emitted graph",
                                {'V', "vanish-consensus"});
    args::Flag _preserve_unmerged_consensus(parser, "bool",
                                            "do not delete original consensus sequences in the merged MAF blocks and in the smoothed graph",
                                            {'N', "preserve-unmerged-consensus"});
    args::ValueFlag<double> _contiguous_path_jaccard(parser, "float",
                                                     "minimum fraction of paths that have to be contiguous for merging MAF blocks and consensus sequences (default: 1.0)",
                                                     {'J', "contiguous-path-jaccard"});

    args::Flag write_block_to_split_fastas(parser, "bool", "write the FASTA sequences for split blocks",
                                           {'A', "write-split-block-fastas"});
    args::Flag write_block_fastas(parser, "bool", "write the FASTA sequences for blocks put into poa",
                                  {'B', "write-poa-block-fastas"});

    args::ValueFlag<std::string> base(parser, "BASE", "use this basename for temporary files during build",
                                      {'b', "base"});
    args::Flag no_prep(parser, "bool",
                       "do not prepare the graph for processing (prep is equivalent to odgi chop followed by odgi sort -p sYgs, and is disabled when taking XG input)",
                       {'n', "no-prep"});
    args::ValueFlag<std::string> _max_block_weight(parser, "N", "maximum seed sequence in block (1k = 1K = 1000, 1m = 1M = 10^6, 1g = 1G = 10^9) [default: 10M]",
                                                {'w', "block-weight-max"});
    args::ValueFlag<std::string> _max_block_jump(parser, "N", "maximum path jump to include in block (1k = 1K = 1000, 1m = 1M = 10^6, 1g = 1G = 10^9) [default: 100]",
                                              {'j', "path-jump-max"});
    args::ValueFlag<std::string> _max_edge_jump(parser, "N", "maximum edge jump before breaking (1k = 1K = 1000, 1m = 1M = 10^6, 1g = 1G = 10^9) [default: 0 / off]",
                                             {'e', "edge-jump-max"});

    args::ValueFlag<uint64_t> _max_merged_groups_in_memory(parser, "N",
                                                           "increasing this value, much more blocks that are not immediately contiguous along the graph will be merged [default: 50]",
                                                           {'G', "max-block-groups-in-memory"});

    // Block split
    args::ValueFlag<double> _block_group_identity(parser, "N",
                                                  "minimum edit-based identity to cluster sequences [default: 0.0]",
                                                  {'I', "block-id-min"});
    args::ValueFlag<double> _block_length_ratio_min(parser, "N",
                                                    "minimum small / large length ratio to cluster in a block [default: 0.0]",
                                                    {'R', "block-ratio-min"});
    args::ValueFlag<std::string> _min_dedup_depth_for_block_splitting(parser, "N",
                                                                   "minimum (deduplicated) block depth for applying the block split (1k = 1K = 1000, 1m = 1M = 10^6, 1g = 1G = 10^9) [default: 0, disabled]",
                                                                   {'d', "min-block-depth-split"});
    args::ValueFlag<std::string> _min_dedup_depth_for_mash_clustering(parser, "N",
                                                                   "minimum (deduplicated) block depth for applying the mash-based clustering (1k = 1K = 1000, 1m = 1M = 10^6, 1g = 1G = 10^9) [default: 12000, 0 to disable it]",
                                                                   {'D', "min-block-depth-mash"});
    args::ValueFlag<std::string> _min_length_mash_based_clustering(parser, "N",
                                                                "minimum sequence length to cluster sequences using mash-distance (1k = 1K = 1000, 1m = 1M = 10^6, 1g = 1G = 10^9) [default: 200, 0 to disable it]",
                                                                {'L', "min-seq-len-mash"});
    args::ValueFlag<double> _block_group_est_identity(parser, "N",
                                                      "minimum mash-based estimated identity to cluster sequences [default: equals to block-id-min]",
                                                      {'E', "block-est-id-max"});

    args::ValueFlag<uint64_t> _kmer_size(parser, "N", "kmer size to compute the mash distance [default: 17]",
                                         {'k', "kmer-size-mash-distance"});

    args::ValueFlag<std::string> _min_copy_length(parser, "N", "minimum repeat length to collapse (1k = 1K = 1000, 1m = 1M = 10^6, 1g = 1G = 10^9) [default: 1000]",
                                               {'c', "copy-length-min"});
    args::ValueFlag<std::string> _max_copy_length(parser, "N",
                                               "maximum repeat length to attempt to detect (1k = 1K = 1000, 1m = 1M = 10^6, 1g = 1G = 10^9) [default: 20K]",
                                               {'W', "copy-length-max"});
    args::ValueFlag<std::string> _target_poa_length(parser, "N", "target length to put into POA, blocks are split when paths go over this length (1k = 1K = 1000, 1m = 1M = 10^6, 1g = 1G = 10^9) [default: 5000]",
                                              {'l', "poa-length-target"});
    args::ValueFlag<std::string> _max_poa_length(parser, "N", "maximum sequence length to put into POA, cut sequences over this length (1k = 1K = 1000, 1m = 1M = 10^6, 1g = 1G = 10^9) [default: 2*poa-length-target = 10k]",
                                              {'q', "poa-length-max"});
    args::ValueFlag<float> _poa_padding_fraction(parser, "N", "flanking sequence length fraction (padding = average sequence length in the block * N) to pad each end of each sequence with during POA, in effect overlapping and trimming the POA problems [default: 0.01]",
                                           {'O', "poa-padding-ratio"});
    args::ValueFlag<std::string> _max_block_depth_for_padding_more(parser, "N",
                                                                   "maximum block depth beyond which a (small) fixed amount of flanking nucleotides is not added (1k = 1K = 1000, 1m = 1M = 10^6, 1g = 1G = 10^9) [default: 1000, 0 to disable it]",
                                                                      {'Y', "max-block-depth-adaptive-poa-padding"});
    args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    args::ValueFlag<uint64_t> num_poa_threads(parser, "N", "use this many POA threads (can be used to reduce memory requirements with large --poa-length-target settings) [default: --threads]", {'T', "poa-threads"});
    args::ValueFlag<std::string> poa_params(parser, "match,mismatch,gap1,ext1(,gap2,ext2)",
                                            "score parameters for partial order alignment, if 4 then gaps are affine, if 6 then gaps are convex [default: 1,4,6,2,26,1]",
                                            {'p', "poa-params"});
    args::ValueFlag<int> _prep_node_chop(parser, "N", "during prep, chop nodes to this length [default: 100]",
                                         {'X', "chop-to"});
    args::ValueFlag<float> _prep_sgd_min_term_updates(parser, "N",
                                                      "path-guided SGD sort quality parameter (N * sum_path_length updates per iteration) for graph prep [default: 1]",
                                                      {'U', "path-sgd-term-updates"});
    args::Flag use_spoa(parser, "use-spoa",
                        "run spoa (in local alignment mode) instead of abPOA (in global alignment mode) for smoothing",
                        {'S', "spoa"});
    args::Flag change_alignment_mode(parser, "change-alignment-mode",
                                     "change the alignment mode of spoa to global, and of abpoa to local",
                                     {'Z', "change-alignment-mode"});
    args::Flag keep_temp(parser, "keep-temp", "keep temporary files", {'K', "keep-temp"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc == 1) {
        std::cout << parser;
        return 1;
    }

    size_t n_threads = num_threads ? args::get(num_threads) : 1;
    omp_set_num_threads(n_threads);
    size_t n_poa_threads = num_poa_threads ? args::get(num_poa_threads) : n_threads;

    std::string smoothed_out_gfa = args::get(smoothed_out);
    std::vector<std::string> consensus_path_names;
    std::vector<smoothxg::consensus_spec_t> consensus_specs;
    bool requires_consensus = !args::get(vanish_consensus);
    bool write_consensus_graph = false;
    std::string consensus_path_prefix = _consensus_path_prefix ? args::get(_consensus_path_prefix) : "Consensus_";

    if (_consensus_spec) {
        consensus_specs = smoothxg::parse_consensus_spec(args::get(_consensus_spec), requires_consensus);
        write_consensus_graph = true;
    }

    if (!_read_consensus_path_names) {
        bool add_consensus = false;
        if (_write_consensus_path_names) {
            add_consensus = true;
        }

        if (requires_consensus) {
            add_consensus = true;
        }

        if (args::get(merge_blocks) && (!write_msa_in_maf_format && !add_consensus)) {
            std::cerr
                    << "[smoothxg::main] error: Please specify -m/--write-msa-in-maf-format and/or keep the consensus "
                       "sequences in the smoothed graph to use the -M/--merge-blocks option." << std::endl;
            return 1;
        }

        if (!args::get(merge_blocks) && (_contiguous_path_jaccard || _preserve_unmerged_consensus)) {
            std::cerr << "[smoothxg::main] error: Please specify -M/--merge-blocks option to use the "
                         "-J/--contiguous-path-jaccard and/or the -N/--preserve-unmerged-consensus option."
                      << std::endl;
            return 1;
        }

        if (args::get(keep_temp) && args::get(no_prep)) {
            std::cerr << "[smoothxg::main] error: Please specify -K/--keep-temp or -n/--no-prep, not both."
                      << std::endl;
            return 1;
        }

        if (!smoothed_out) {
            std::cerr << "[smoothxg::main] error: Please specify an output file with -o/--smoothed-out." << std::endl;
            return 1;
        }


        const double contiguous_path_jaccard = _contiguous_path_jaccard ? min(args::get(_contiguous_path_jaccard), 1.0) : 1.0;

        const uint64_t max_block_weight = _max_block_weight ? (uint64_t)smoothxg::handy_parameter(args::get(_max_block_weight), 10000000) : 10000000;
        const uint64_t max_block_jump = _max_block_jump ? (uint64_t)smoothxg::handy_parameter(args::get(_max_block_jump), 100) : 100;
        const uint64_t max_edge_jump = _max_edge_jump ? (uint64_t)smoothxg::handy_parameter(args::get(_max_edge_jump), 0) : 0;
        const uint64_t min_copy_length = _min_copy_length ? (uint64_t)smoothxg::handy_parameter(args::get(_min_copy_length), 1000) : 1000;
        const uint64_t max_copy_length = _max_copy_length ? (uint64_t)smoothxg::handy_parameter(args::get(_max_copy_length), 20000) : 20000;
        const uint64_t target_poa_length = _target_poa_length ? (uint64_t)smoothxg::handy_parameter(args::get(_target_poa_length), 5000) : 5000;
        const uint64_t max_poa_length = _max_poa_length ? (uint64_t)smoothxg::handy_parameter(args::get(_max_poa_length), 2 * target_poa_length) : 2 * target_poa_length;
        const float poa_padding_fraction = _poa_padding_fraction ? args::get(_poa_padding_fraction) : 0.01;
        const uint64_t max_block_depth_for_padding_more = _max_block_depth_for_padding_more ?
                (uint64_t)smoothxg::handy_parameter(args::get(_max_block_depth_for_padding_more), 1000) : 1000;

        const uint64_t max_merged_groups_in_memory = _max_merged_groups_in_memory ? args::get(_max_merged_groups_in_memory)
                                                                            : 50;

        // Block split
        const double block_length_ratio_min = _block_length_ratio_min ? args::get(_block_length_ratio_min) : 0.0;
        const uint64_t min_length_mash_based_clustering = _min_length_mash_based_clustering ?
                (uint64_t)smoothxg::handy_parameter(args::get(_min_length_mash_based_clustering), 200) : 200;
        const uint64_t kmer_size = _kmer_size ? args::get(_kmer_size) : 17;
        if (min_length_mash_based_clustering != 0 && min_length_mash_based_clustering < kmer_size) {
            std::cerr
                    << "[smoothxg::main] error: the minimum sequences length to cluster sequences using mash-distance "
                       "has to be greater than or equal to the kmer size."
                    << std::endl;
            return 1;
        }

        const uint64_t min_dedup_depth_for_block_splitting = _min_dedup_depth_for_block_splitting ?
                (uint64_t)smoothxg::handy_parameter(args::get(_min_dedup_depth_for_block_splitting), 0) : 0;

        const double block_group_identity = _block_group_identity ? args::get(_block_group_identity) : 0.0;
        const double block_group_est_identity = _block_group_est_identity ? args::get(_block_group_est_identity)
                                                                    : block_group_identity;
        const uint64_t min_dedup_depth_for_mash_clustering = _min_dedup_depth_for_mash_clustering ?
                (uint64_t)smoothxg::handy_parameter(args::get(_min_dedup_depth_for_mash_clustering), 12000) : 12000;

        int poa_m = 1;
        int poa_n = 4;
        int poa_g = 6;
        int poa_e = 2;
        int poa_q = 26;
        int poa_c = 1;

        if (!args::get(poa_params).empty()) {
            if (args::get(poa_params).find(',') == std::string::npos) {
                std::cerr
                        << "[smoothxg::main] error: either 4 or 6 POA scoring parameters must be given to -p --poa-params"
                        << std::endl;
                return 1;
            }
            std::vector<std::string> params_str = smoothxg::split(args::get(poa_params), ',');
            std::vector<int> params(params_str.size());
            std::transform(params_str.begin(), params_str.end(), params.begin(),
                           [](const std::string &s) { return std::stoi(s); });
            if (params.size() == 6) {
                poa_m = params[0];
                poa_n = params[1];
                poa_g = params[2];
                poa_e = params[3];
                poa_q = params[4];
                poa_c = params[5];
            } else if (params.size() == 4) {
                poa_m = params[0];
                poa_n = params[1];
                poa_g = params[2];
                poa_e = params[3];
                if (args::get(use_spoa)) {
                    poa_q = poa_g;
                    poa_c = poa_e;
                } else {
                    poa_q = 0;
                    poa_c = 0;
                }
            } else {
                std::cerr
                        << "[smoothxg::main] error: either 4 or 6 POA scoring parameters must be given to -p --poa-params"
                        << std::endl;
                return 1;
            }

        }

        const bool order_paths_from_longest = true; //args::get(use_spoa);
        const float term_updates = (_prep_sgd_min_term_updates ? args::get(_prep_sgd_min_term_updates) : 1);
        const int node_chop = (_prep_node_chop ? args::get(_prep_node_chop) : 100);

        std::cerr << "[smoothxg::main] loading graph" << std::endl;
        auto graph = std::make_unique<XG>();
        if (!args::get(xg_in).empty()) {
            std::ifstream in(args::get(xg_in));
            graph->deserialize(in);
        } else if (!args::get(gfa_in).empty()) {
            // prep the graph by default
            std::string gfa_in_name;
            if (!args::get(no_prep)) {
                if (args::get(base).empty()) {
                    gfa_in_name = args::get(gfa_in) + ".prep.gfa";
                } else {
                    gfa_in_name = args::get(base) + '/' + args::get(gfa_in) + ".prep.gfa";
                }
                std::cerr << "[smoothxg::main] prepping graph for smoothing" << std::endl;
                smoothxg::prep(args::get(gfa_in), gfa_in_name, node_chop, term_updates, true, args::get(base), n_threads);
            } else {
                gfa_in_name = args::get(gfa_in);
            }
            std::cerr << "[smoothxg::main] building xg index" << std::endl;
            graph->from_gfa(gfa_in_name, false, args::get(base));
            if (!args::get(keep_temp) && !args::get(no_prep)) {
                std::remove(gfa_in_name.c_str());
            }
        }

        auto *blockset = new smoothxg::blockset_t();
        smoothxg::smoothable_blocks(*graph,
                                    *blockset,
                                    max_block_weight,
                                    target_poa_length,
                                    max_block_jump,
                                    max_edge_jump,
                                    order_paths_from_longest,
                                    num_threads);

        const uint64_t min_autocorr_z = 5;
        const uint64_t autocorr_stride = 50;

        smoothxg::break_blocks(*graph,
                               blockset,
                               block_length_ratio_min,
                               min_length_mash_based_clustering,
                               block_group_identity,
                               block_group_est_identity,
                               kmer_size,
                               min_dedup_depth_for_block_splitting,
                               min_dedup_depth_for_mash_clustering,
                               max_poa_length,
                               min_copy_length,
                               max_copy_length,
                               min_autocorr_z,
                               autocorr_stride,
                               order_paths_from_longest,
                               true,
                               n_threads,
                               args::get(write_block_to_split_fastas));

        // build the path_step_rank_ranges -> index_in_blocks_vector
        // flat_hash_map using SKA: KEY: path_name, VALUE: sorted interval_tree using cgranges https://github.com/lh3/cgranges:
        // we collect path_step_rank_ranges and the identifier of an interval is the index of a block in the blocks vector
        //ska::flat_hash_map<std::string, IITree<uint64_t , uint64_t>> happy_tree_friends = smoothxg::generate_path_nuc_range_block_index(blocks, graph);

        const bool local_alignment = args::get(use_spoa) ^args::get(change_alignment_mode);

        std::string maf_header;
        if (write_msa_in_maf_format) {
            basic_string<char> filename;
            if (!args::get(xg_in).empty()) {
                size_t found = args::get(xg_in).find_last_of("/\\");
                filename = (args::get(xg_in).substr(found + 1));
            } else if (!args::get(gfa_in).empty()) {
                size_t found = args::get(gfa_in).find_last_of("/\\");
                filename = (args::get(gfa_in).substr(found + 1));
            }

            maf_header += "##maf version=1\n";
            maf_header += "# smoothxg\n";
            maf_header += "# input=" + filename + " sequences=" + std::to_string(graph->get_path_count()) + "\n";

            // Merge mode
            maf_header += "# merge_blocks=";
            maf_header += (args::get(merge_blocks) ? "true" : "false");
            maf_header += " contiguous_path_jaccard=" + std::to_string(contiguous_path_jaccard) + "\n";

            // POA
            maf_header += "# POA=";
            maf_header += (args::get(use_spoa) ? "SPOA" : "abPOA");
            maf_header += " alignment_mode=";
            maf_header += (local_alignment ? "local" : "global");
            maf_header += " order_paths=from_";
            maf_header += (order_paths_from_longest ? "longest" : "shortest");
            maf_header += "\n";

            // create_blocks
            maf_header += "# max_block_weight=" + std::to_string(max_block_weight) +
                          " max_block_jump=" + std::to_string(max_block_jump) +
                          " max_edge_jump=" + std::to_string(max_edge_jump) + "\n";

            // break_blocks
            maf_header += "# max_poa_length=" + std::to_string(max_poa_length) +
                          " min_copy_length=" + std::to_string(min_copy_length) +
                          " max_copy_length=" + std::to_string(max_copy_length) +
                          " min_autocorr_z=" + std::to_string(min_autocorr_z) +
                          " autocorr_stride=" + std::to_string(autocorr_stride) + "\n";

            // split_blocks
            maf_header += "# block_group_identity=" + std::to_string(block_group_identity) +
                          " block_group_estimated_identity=" + std::to_string(block_group_est_identity) +
                          " min_length_mash_based_clustering=" + std::to_string(min_length_mash_based_clustering) +
                          " min_dedup_depth_for_mash_clustering=" +
                          std::to_string(min_dedup_depth_for_mash_clustering) +
                          " kmer_size=" + std::to_string(_kmer_size) + "\n";
        }

        {
            auto smoothed = smoothxg::smooth_and_lace(*graph,
                                                      blockset,
                                                      poa_m,
                                                      poa_n,
                                                      poa_g,
                                                      poa_e,
                                                      poa_q,
                                                      poa_c,
                                                      poa_padding_fraction,
                                                      max_block_depth_for_padding_more,
                                                      local_alignment,
                                                      n_threads,
                                                      n_poa_threads,
                                                      args::get(write_msa_in_maf_format), maf_header,
                                                      args::get(merge_blocks), args::get(_preserve_unmerged_consensus),
                                                      contiguous_path_jaccard,
                                                      !args::get(use_spoa),
                                                      add_consensus ? consensus_path_prefix : "",
                                                      consensus_path_names,
                                                      args::get(write_block_fastas),
                                                      max_merged_groups_in_memory);

            std::cerr << "[smoothxg::main] unchopping smoothed graph" << std::endl;
            odgi::algorithms::unchop(*smoothed, n_threads, true);

            uint64_t smoothed_nodes = 0;
            uint64_t smoothed_length = 0;
            smoothed->for_each_handle(
                    [&](const handle_t &h) {
                        ++smoothed_nodes;
                        smoothed_length += smoothed->get_length(h);
                    });
            std::cerr << "[smoothxg::main] smoothed graph length " << smoothed_length << "bp " << "in "
                      << smoothed_nodes << " nodes" << std::endl;

            std::cerr << "[smoothxg::main] writing smoothed graph to " << smoothed_out_gfa << std::endl;
            ofstream out(smoothed_out_gfa.c_str());
            smoothed->to_gfa(out);
            out.close();
            delete smoothed;
        }

        delete blockset;

        // do we need to write the consensus path names?
        if (_write_consensus_path_names) {
            std::string write_consensus_path_names = args::get(_write_consensus_path_names);
            std::ofstream consensus_path_names_out(write_consensus_path_names);
            for (auto &consensus_path_name : consensus_path_names) {
                consensus_path_names_out << consensus_path_name << std::endl;
            }
            consensus_path_names_out.close();
        }

        // end !_read_consenus_path_names
    } else {
        if (!_smoothed_in_gfa) {
            std::cerr << "[smoothxg::main] error: Please only use the -H/--read-consensus-path-names parameter"
                         " together with the -F/--smoothed-in option." << std::endl;
            return 1;
        }
    }

    // do we need to build the consensus graph?
    if (write_consensus_graph) {
        // check if we have a reference path list
        /*
        if (_ref_paths) {
            if (_only_ref_paths) {
                consensus_path_names.clear();
            }
        }
        // get the base name
        std::string consensus_base = args::get(write_consensus_graph);
        std::vector<uint64_t> jump_maxes;
        if (_consensus_jump_max) {
            for (auto& s : smoothxg::split(args::get(_consensus_jump_max),',')) {
                jump_maxes.push_back(std::stoi(s));
            }
        } else {
            jump_maxes.push_back(100);
        }
        */
        //uint64_t jump_limit = (_consensus_jump_limit ? args::get(_consensus_jump_limit) : 1e6);
        std::cerr << "[smoothxg::main] building xg index from smoothed graph" << std::endl;
        XG smoothed_xg;
        if (_read_consensus_path_names) {
            std::string smoothed_in_gfa = args::get(_smoothed_in_gfa);
            smoothed_xg.from_gfa(smoothed_in_gfa, false,
                                 args::get(base).empty() ? smoothed_in_gfa : args::get(base));
            std::ifstream file(args::get(_read_consensus_path_names));
            std::string path_name;
            while (std::getline(file, path_name)) {
                consensus_path_names.push_back(path_name);
            }
        } else {
            smoothed_xg.from_gfa(smoothed_out_gfa, false,
                                 args::get(base).empty() ? smoothed_out_gfa : args::get(base));
        }
        for (auto &spec : consensus_specs) {
            //for (auto jump_max : jump_maxes) {
            std::vector<std::string> consensus_paths_to_use;
            if (spec.ref_file.size()) {
                ifstream ref_paths(spec.ref_file.c_str());
                std::string line;
                while (std::getline(ref_paths, line)) {
                    consensus_paths_to_use.push_back(line);
                }
            }
            if (spec.keep_consensus_paths) {
                consensus_paths_to_use.insert(consensus_paths_to_use.begin(),
                                              consensus_path_names.begin(),
                                              consensus_path_names.end());
            }
            auto outname = displayname(spec) + ".gfa";
            std::cerr << "[smoothxg::create_consensus_graph] deriving consensus graph " << outname << std::endl;
            odgi::graph_t *consensus_graph = smoothxg::create_consensus_graph(smoothed_xg,
                                                                              consensus_paths_to_use,
                                                                              spec.min_allele_len,
                                                                              spec.max_allele_len,
                                                                              spec.min_consensus_path_cov,
                                                                              n_threads,
                                                                              outname);
            ofstream o(outname);
            consensus_graph->to_gfa(o);
            o.close();
            delete consensus_graph;
        }
    }

    return 0;
}
