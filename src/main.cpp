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
#include <chrono>

using namespace std;
using namespace xg;

/*
std::string print_time(const double &_seconds) {
    int days = 0, hours = 0, minutes = 0, seconds = 0;
    distribute_seconds(days, hours, minutes, seconds, _seconds);
    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(2) << days << ":"
           << std::setfill('0') << std::setw(2) << hours << ":"
           << std::setfill('0') << std::setw(2) << minutes << ":"
           << std::setfill('0') << std::setw(2) << seconds;
    return buffer.str();
}

void distribute_seconds(int &days, int &hours, int &minutes, int &seconds, const double &input_seconds) {
    const int cseconds_in_day = 86400;
    const int cseconds_in_hour = 3600;
    const int cseconds_in_minute = 60;
    const int cseconds = 1;
    days = std::floor(input_seconds / cseconds_in_day);
    hours = std::floor(((int) input_seconds % cseconds_in_day) / cseconds_in_hour);
    minutes = std::floor((((int) input_seconds % cseconds_in_day) % cseconds_in_hour) / cseconds_in_minute);
    seconds = ((((int) input_seconds % cseconds_in_day) % cseconds_in_hour) % cseconds_in_minute) /
              cseconds;
}
 */


int main(int argc, char** argv) {
    args::ArgumentParser parser("smoothxg: collinear block finder and graph consensus generator");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> gfa_in(parser, "FILE", "index the graph in this GFA file", {'g', "gfa-in"});
    args::ValueFlag<std::string> xg_in(parser, "FILE", "read the xg index from this file", {'i', "in"});
    args::ValueFlag<std::string> smoothed_out(parser, "FILE", "write GFA to this file (not /dev/stdout if consensus graph is made)", {'o', "smoothed-out"});
    args::ValueFlag<std::string> _smoothed_in_gfa(parser, "FILE", "read GFA from this file", {'F', "smoothed-in"});
    args::ValueFlag<std::string> write_msa_in_maf_format(parser, "FILE","write the multiple sequence alignments (MSAs) in MAF format in this file",{'m', "write-msa-in-maf-format"});
    args::Flag add_consensus(parser, "bool", "include consensus sequence in the smoothed graph", {'a', "add-consensus"});
    args::ValueFlag<std::string> _write_consensus_path_names(parser, "FILE", "write the consensus path names to this file", {'f', "write-consensus-path-names"});
    args::ValueFlag<std::string> _read_consensus_path_names(parser, "FILE", "read the consensus path names from this file", {'D', "read-consensus-path-names"});
    args::ValueFlag<std::string> write_consensus_graph(parser, "BASENAME", "write the consensus graph to BASENAME.cons_[jump_max].gfa", {'s', "write-consensus-graph"});
    args::ValueFlag<std::string> _consensus_jump_max(parser, "jump_max[,jump_max]*", "preserve all divergences from the consensus paths greater than this length, with multiples allowed [default: 100]", {'C', "consensus-jump-max"});

    // Merge blocks (for merging MAF blocks and consensus sequences)
    args::Flag merge_blocks(parser, "bool", "merge contiguous MAF blocks in the MAF output and consensus sequences in the smoothed graph",{'M', "merge-blocks"});
    args::Flag _preserve_unmerged_consensus(parser, "bool", "do not delete original consensus sequences in the merged MAF blocks and in the smoothed graph",{'N', "preserve-unmerged-consensus"});
    args::ValueFlag<double> _contiguous_path_jaccard(parser, "float","minimum fraction of paths that have to be contiguous for merging MAF blocks and consensus sequences (default: 1.0)",{'J', "contiguous-path-jaccard"});

    args::Flag write_block_to_split_fastas(parser, "bool", "write the FASTA sequences for split blocks",{'A', "write-split-block-fastas"});
    args::Flag write_block_fastas(parser, "bool", "write the FASTA sequences for blocks put into poa",{'B', "write-poa-block-fastas"});

    args::ValueFlag<std::string> base(parser, "BASE", "use this basename for temporary files during build", {'b', "base"});
    args::Flag no_prep(parser, "bool", "do not prepare the graph for processing (prep is equivalent to odgi chop followed by odgi sort -p sYgs, and is disabled when taking XG input)", {'n', "no-prep"});
    args::ValueFlag<uint64_t> _max_block_weight(parser, "N", "maximum seed sequence in block [default: 10000]", {'w', "block-weight-max"});
    args::ValueFlag<uint64_t> _max_block_jump(parser, "N", "maximum path jump to include in block [default: 5000]", {'j', "path-jump-max"});
    args::ValueFlag<uint64_t> _min_subpath(parser, "N", "minimum length of a subpath to include in partial order alignment [default: 0 / no filter]", {'k', "subpath-min"});
    args::ValueFlag<uint64_t> _max_edge_jump(parser, "N", "maximum edge jump before breaking [default: 5000]", {'e', "edge-jump-max"});
    //args::ValueFlag<double> _min_segment_ratio(parser, "N", "split out segments in a block that are less than this fraction of the length of the longest path range in the block [default: 0.1]", {'R', "min-segment-ratio"});
    args::ValueFlag<double> _block_group_identity(parser, "N", "minimum edit-based identity to group sequences in POA [default: 0.5]", {'I', "block-id-min"});
    args::ValueFlag<uint64_t> _min_copy_length(parser, "N", "minimum repeat length to collapse [default: 1000]", {'c', "copy-length-min"});
    args::ValueFlag<uint64_t> _max_copy_length(parser, "N", "maximum repeat length to attempt to detect [default: 20000]", {'W', "copy-length-max"});
    args::ValueFlag<uint64_t> _max_poa_length(parser, "N", "maximum sequence length to put into poa [default: 10000]", {'l', "poa-length-max"});
    args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    args::ValueFlag<std::string> poa_params(parser, "match,mismatch,gap1,ext1(,gap2,ext2)", "score parameters for partial order alignment, if 4 then gaps are affine, if 6 then gaps are convex [default: 2,4,4,2,24,1]", {'p', "poa-params"});
    args::ValueFlag<int> _prep_node_chop(parser, "N", "during prep, chop nodes to this length [default: 100]", {'X', "chop-to"});
    args::ValueFlag<float> _prep_sgd_min_term_updates(parser, "N", "path-guided SGD sort quality parameter (N * sum_path_length updates per iteration) for graph prep [default: 1]", {'U', "path-sgd-term-updates"});
    args::Flag use_spoa(parser, "use-spoa", "run spoa (in local alignment mode) instead of abPOA (in global alignment mode) for smoothing", {'S', "spoa"});
    args::Flag change_alignment_mode(parser, "change-alignment-mode", "change the alignment mode of spoa to global, the local alignment mode of abpoa is currently not supported", {'Z', "change-alignment-mode"});
    args::Flag no_toposort(parser, "prep-no-toposort", "do not apply topological sorting in the prep sort pipeline", {'T', "prep-no-toposort"});
    args::Flag validate(parser, "validate", "validate construction", {'V', "validate"});
    args::Flag keep_temp(parser, "keep-temp", "keep temporary files", {'K', "keep-temp"});
    args::Flag debug(parser, "debug", "enable debugging", {'d', "debug"});
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
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    size_t n_threads = num_threads ? args::get(num_threads) : 1;
    omp_set_num_threads(n_threads);

    std::string smoothed_out_gfa = args::get(smoothed_out);
    std::vector<std::string> consensus_path_names;

    if (!_read_consensus_path_names) {

    if (args::get(merge_blocks) && (!write_msa_in_maf_format && !args::get(add_consensus))) {
        std::cerr << "[smoothxg::main] error: Please specify -m/--write-msa-in-maf-format and/or -a/--add-consensus "
                     "to use the -M/--merge-blocks option." << std::endl;
        return 1;
    }

    if (!args::get(add_consensus) && write_consensus_graph) {
        std::cerr << "[smoothxg::main] error: Please only use the -s/--write-consensus-graph parameter together with"
                   "the -a/--add-consensus option." << std::endl;
        return 1;
    }

    if (_min_subpath && write_consensus_graph) {
        std::cerr << "[smoothxg::main] error: Please only use the -s/--write-consensus-graph parameter without"
                   "the -k/--subpath option." << std::endl;
        return 1;
    }

    if (!args::get(merge_blocks) && (_contiguous_path_jaccard || _preserve_unmerged_consensus)) {
        std::cerr << "[smoothxg::main] error: Please specify -M/--merge-blocks option to use the "
                     "-J/--contiguous-path-jaccard and/or the -N/--preserve-unmerged-consensus option." << std::endl;
        return 1;
    }

    if (args::get(keep_temp) && args::get(no_prep)) {
        std::cerr << "[smoothxg::main] error: Please specify -K/--keep-temp or -n/--no-prep, not both." << std::endl;
        return 1;
    }

    if (!smoothed_out) {
        std::cerr << "[smoothxg::main] error: Please specify an output file with -o/--smoothed-out." << std::endl;
        return 1;
    }

    if (_write_consensus_path_names && !add_consensus) {
        std::cerr << "[smoothxg::main] error: Please use -f/--write-consensus-path-names only with -a/--add-consensus." << std::endl;
        return 1;
    }
    double contiguous_path_jaccard = _contiguous_path_jaccard ? min(args::get(_contiguous_path_jaccard), 1.0) : 1.0;

    uint64_t max_block_weight = _max_block_weight ? args::get(_max_block_weight) : 10000;
    uint64_t max_block_jump = _max_block_jump ? args::get(_max_block_jump) : 5000;
    uint64_t min_subpath = _min_subpath ? args::get(_min_subpath) : 0;
    uint64_t max_edge_jump = _max_edge_jump ? args::get(_max_edge_jump) : 5000;
    uint64_t min_copy_length = _min_copy_length ? args::get(_min_copy_length) : 1000;
    uint64_t max_copy_length = _max_copy_length ? args::get(_max_copy_length) : 20000;
    uint64_t max_poa_length = _max_poa_length ? args::get(_max_poa_length) : 10000;
    double block_group_identity =  _block_group_identity ? args::get(_block_group_identity) : 0.5;

    if (!args::get(use_spoa) && args::get(change_alignment_mode)) {
        std::cerr
                << "[smoothxg::main] error: Currently, the local alignment mode of abpoa is not supported. As default "
                << "abpoa is ran in global mode. You can select spoa in local alignment mode via -S, --spoa. To run spoa in "
                   "global mode, please additionally specify -Z, --change-alignment-mode."
                << std::endl;
        return 1;
    }

    int poa_m = 2;
    int poa_n = 4;
    int poa_g = 4;
    int poa_e = 2;
    int poa_q = 24;
    int poa_c = 1;

    if (!args::get(poa_params).empty()) {
        if (args::get(poa_params).find(',') == std::string::npos) {
            std::cerr << "[smoothxg::main] error: either 4 or 6 POA scoring parameters must be given to -p --poa-params" << std::endl;
            return 1;
        }
        std::vector<std::string> params_str = smoothxg::split(args::get(poa_params),',');
        std::vector<int> params(params_str.size());
        std::transform(params_str.begin(), params_str.end(), params.begin(),
                       [](const std::string& s) { return std::stoi(s); });
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
            std::cerr << "[smoothxg::main] error: either 4 or 6 POA scoring parameters must be given to -p --poa-params" << std::endl;
            return 1;
        }

    }

    bool order_paths_from_longest = args::get(use_spoa);
    float term_updates = (_prep_sgd_min_term_updates ? args::get(_prep_sgd_min_term_updates) : 1);
    float node_chop = (_prep_node_chop ? args::get(_prep_node_chop) : 100);

    std::cerr << "[smoothxg::main] loading graph" << std::endl;
    auto graph = std::make_unique<XG>();
    if (!args::get(xg_in).empty()) {
        std::ifstream in(args::get(xg_in));
        graph->deserialize(in);
    } else if (!args::get(gfa_in).empty()) {
        // prep the graph by default
        std::string gfa_in_name;
        if (!args::get(no_prep)) {
            if (args::get(base).empty()){
                gfa_in_name = args::get(gfa_in) + ".prep.gfa";
            }else{
                gfa_in_name = args::get(base) + '/' + args::get(gfa_in) + ".prep.gfa";
            }
            std::cerr << "[smoothxg::main] prepping graph for smoothing" << std::endl;
            smoothxg::prep(args::get(gfa_in), gfa_in_name, node_chop, term_updates, !args::get(no_toposort), n_threads);
        } else {
            gfa_in_name = args::get(gfa_in);
        }
        std::cerr << "[smoothxg::main] building xg index" << std::endl;
        graph->from_gfa(gfa_in_name, args::get(validate),
                       args::get(base).empty() ? gfa_in_name : args::get(base));
        if (!args::get(keep_temp) && !args::get(no_prep)) {
            std::remove(gfa_in_name.c_str());
        }
    }

    auto* blockset = new smoothxg::blockset_t("blocks");
    smoothxg::smoothable_blocks(*graph,
                                *blockset,
                                max_block_weight,
                                max_block_jump,
                                min_subpath,
                                max_edge_jump,
                                order_paths_from_longest,
                                num_threads);

    uint64_t min_autocorr_z = 5;
    uint64_t autocorr_stride = 50;

    smoothxg::break_blocks(*graph,
                           blockset,
                           block_group_identity,
                           max_poa_length,
                           min_copy_length,
                           max_copy_length,
                           min_autocorr_z,
                           autocorr_stride,
                           order_paths_from_longest,
                           true,
                           n_threads,
                           write_consensus_graph,
                           args::get(write_block_to_split_fastas));

    // build the path_step_rank_ranges -> index_in_blocks_vector
    // flat_hash_map using SKA: KEY: path_name, VALUE: sorted interval_tree using cgranges https://github.com/lh3/cgranges:
    // we collect path_step_rank_ranges and the identifier of an interval is the index of a block in the blocks vector
    //ska::flat_hash_map<std::string, IITree<uint64_t , uint64_t>> happy_tree_friends = smoothxg::generate_path_nuc_range_block_index(blocks, graph);

    bool local_alignment = args::get(use_spoa) ^ args::get(change_alignment_mode);

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
                " min_subpath=" + std::to_string(min_subpath) +
                " max_edge_jump=" + std::to_string(max_edge_jump) + "\n";

        // break_blocks
        maf_header += "# max_poa_length=" + std::to_string(max_poa_length) +
                " min_copy_length=" + std::to_string(min_copy_length) +
                " max_copy_length=" + std::to_string(max_copy_length) +
                " min_autocorr_z=" + std::to_string(min_autocorr_z) +
                " autocorr_stride=" + std::to_string(autocorr_stride) +
                " block_group_identity=" + std::to_string(block_group_identity) + "\n";
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
                                                  local_alignment,
                                                  n_threads,
                                                  args::get(write_msa_in_maf_format), maf_header,
                                                  args::get(merge_blocks), args::get(_preserve_unmerged_consensus), contiguous_path_jaccard,
                                                  !args::get(use_spoa),
                                                  args::get(add_consensus) ? "Consensus_" : "",
                                                  consensus_path_names,
                                                  args::get(write_block_fastas));

        uint64_t smoothed_nodes = 0;
        uint64_t smoothed_length = 0;
        smoothed->for_each_handle(
            [&](const handle_t& h) {
                ++smoothed_nodes;
                smoothed_length += smoothed->get_length(h);
            });
        std::cerr << "[smoothxg::main] smoothed graph length " << smoothed_length << "bp " << "in " << smoothed_nodes << " nodes" << std::endl;

        std::cerr << "[smoothxg::main] writing smoothed graph to " << smoothed_out_gfa << std::endl;
        ofstream out(smoothed_out_gfa.c_str());
        smoothed->to_gfa(out);
        out.close();
        delete smoothed;
        delete blockset;
    }

    // do we need to write the consensus path names?
    if (_write_consensus_path_names) {
        std::string write_consensus_path_names = args::get(_write_consensus_path_names);
        std::ofstream consensus_path_names_out(write_consensus_path_names);
        for (auto& consensus_path_name : consensus_path_names) {
            consensus_path_names_out << consensus_path_name << std::endl;
        }
        consensus_path_names_out.close();
    }

    // end !_read_consenus_path_names
    } else {
        if (!_smoothed_in_gfa) {
            std::cerr << "[smoothxg::main] error: Please only use the -D/--read-consensus-path-names parameter"
                         " together with the -F/--smoothed-in option." << std::endl;
        return 1;
        }
    }

    // do we need to build the consensus graph?
    if (write_consensus_graph) {
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
        std::cerr << "[smoothxg::main] building xg index from smoothed graph" << std::endl;
        XG smoothed_xg;
        if (_read_consensus_path_names) {
            std::string smoothed_in_gfa = args::get(_smoothed_in_gfa);
            smoothed_xg.from_gfa(smoothed_in_gfa, args::get(validate),
                                 args::get(base).empty() ? smoothed_in_gfa : args::get(base));
            std::ifstream file(args::get(_read_consensus_path_names));
            std::string path_name;
            while (std::getline(file, path_name)) {
                consensus_path_names.push_back(path_name);
            }
        } else {
            smoothed_xg.from_gfa(smoothed_out_gfa, args::get(validate),
                                 args::get(base).empty() ? smoothed_out_gfa : args::get(base));
        }
        for (auto jump_max : jump_maxes) {
            odgi::graph_t* consensus_graph = smoothxg::create_consensus_graph(
                smoothed_xg, consensus_path_names, jump_max, n_threads,
                args::get(base).empty() ? args::get(write_consensus_graph) : args::get(base));
            ofstream o(consensus_base + "@" + std::to_string(jump_max) + ".gfa");
            consensus_graph->to_gfa(o);
            o.close();
            delete consensus_graph;
        }
    }

    return 0;
}
