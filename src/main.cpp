/** \file smoothxg
 *
 * smooth a graph
 */


#include <omp.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "args.hxx"
#include "sdsl/bit_vectors.hpp"
#include "chain.hpp"
#include "blocks.hpp"
#include "smooth.hpp"
#include "xg.hpp"
#include "prep.hpp"
#include "odgi/odgi.hpp"

using namespace std;
using namespace xg;

int main(int argc, char** argv) {
    args::ArgumentParser parser("smoothxg: collinear block finder and graph consensus generator");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> gfa_in(parser, "FILE", "index the graph in this GFA file", {'g', "gfa-in"});
    //args::ValueFlag<std::string> xg_out(parser, "FILE", "write the resulting xg index to this file", {'o', "out"});
    args::ValueFlag<std::string> xg_in(parser, "FILE", "read the xg index from this file", {'i', "in"});
    args::ValueFlag<std::string> base(parser, "BASE", "use this basename for temporary files during build", {'b', "base"});
    args::Flag no_prep(parser, "bool", "do not prepare the graph for processing (prep is equivalent to odgi chop followed by odgi sort -p sYgs, and is disabled when taking XG input)", {'n', "no-prep"});
    args::Flag add_consensus(parser, "bool", "include consensus sequence in graph", {'a', "add-consensus"});
    args::ValueFlag<uint64_t> _max_block_weight(parser, "N", "maximum seed sequence in block [default: 10000]", {'w', "max-block-weight"});
    args::ValueFlag<uint64_t> _max_block_jump(parser, "N", "maximum path jump to include in block [default: 1000]", {'j', "max-path-jump"});
    args::ValueFlag<uint64_t> _min_subpath(parser, "N", "minimum length of a subpath to include in partial order alignment [default: 0]", {'k', "min-subpath"});
    args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    args::ValueFlag<int> _poa_m(parser, "N", "spoa score for matching bases [default: 5]", {'M', "poa-match"});
    args::ValueFlag<int> _poa_n(parser, "N", "spoa score for mismatching bases [default: -4]", {'N', "poa-mismatch"});
    args::ValueFlag<int> _poa_g(parser, "N", "spoa gap opening penalty (must be negative) [default: -8]", {'G', "poa-gap-open"});
    args::ValueFlag<int> _poa_e(parser, "N", "spoa gap extension penalty (must be negative) [default: -6]", {'E', "poa-gap-extend"});
    args::ValueFlag<int> _poa_q(parser, "N", "spoa gap opening penalty of the second affine function (must be negative) [default: -10]", {'Q', "poa-2nd-gap-open"});
    args::ValueFlag<int> _poa_c(parser, "N", "spoa gap extension penalty of the second affine function (must be negative) [default: -4]", {'C', "poa-2nd-gap-extend"});
    args::ValueFlag<int> _prep_node_chop(parser, "N", "during prep, chop nodes to this length [default: 32]", {'X', "chop-to"});
    args::ValueFlag<float> _prep_sgd_min_term_updates(parser, "N", "path-guided SGD sort quality parameter (N*graph_size updates per iteration) for graph prep [default: 10]", {'U', "path-sgd-term-updates"});
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

    size_t n_threads = args::get(num_threads);
    if (n_threads) {
        omp_set_num_threads(args::get(num_threads));
    } else {
        omp_set_num_threads(1);
    }

    XG graph;
    if (!args::get(xg_in).empty()) {
        std::ifstream in(args::get(xg_in));
        graph.deserialize(in);
    } else if (!args::get(gfa_in).empty()) {
        // prep the graph by default
        std::string gfa_in_name;
        if (!args::get(no_prep)) {
            gfa_in_name = args::get(gfa_in) + ".prep.gfa";
            float term_updates = (_prep_sgd_min_term_updates ? args::get(_prep_sgd_min_term_updates) : 10);
            float node_chop = (_prep_node_chop ? args::get(_prep_node_chop) : 32);
            smoothxg::prep(args::get(gfa_in), gfa_in_name, node_chop, term_updates);
        } else {
            gfa_in_name = args::get(gfa_in);
        }
        graph.from_gfa(gfa_in_name, args::get(validate),
                       args::get(base).empty() ? gfa_in_name : args::get(base));
        if (!args::get(keep_temp) && !args::get(no_prep)) {
            std::remove(gfa_in_name.c_str());
        }
    }

    uint64_t max_block_weight = _max_block_weight ? args::get(_max_block_weight) : 10000;
    uint64_t max_block_jump = _max_block_jump ? args::get(_max_block_jump) : 1000;
    uint64_t min_subpath = _min_subpath ? args::get(_min_subpath) : 0;

    std::int8_t poa_m = 5;
    std::int8_t poa_n = -4;
    std::int8_t poa_g = -8;
    std::int8_t poa_e = -6;
    std::int8_t poa_q = -10;
    std::int8_t poa_c = -4;

    if (_poa_m) poa_m = args::get(_poa_m);
    if (_poa_n) poa_n = args::get(_poa_n);
    if (_poa_g) poa_g = args::get(_poa_g);
    if (_poa_e) poa_e = args::get(_poa_e);
    if (_poa_q) poa_q = args::get(_poa_q);
    if (_poa_c) poa_c = args::get(_poa_c);

    auto blocks = smoothxg::smoothable_blocks(graph,
                                              max_block_weight,
                                              max_block_jump,
                                              min_subpath);

    auto smoothed = smoothxg::smooth_and_lace(graph,
                                              blocks,
                                              poa_m,
                                              poa_n,
                                              poa_g,
                                              poa_e,
                                              poa_q,
                                              poa_c,
                                              args::get(add_consensus) ? "Consensus_" : "");

    smoothed.to_gfa(std::cout);
    
    /*


    uint64_t block_id = 0;
    for (auto& block : blocks) {
        std::cout << "block" << block_id++ << "\t"
                  << block.total_path_length << "\t"
                  << block.max_path_length << "\t"
                  << graph.get_id(block.handles.front())
                  << "-" << graph.get_id(block.handles.back()) << "\t"
                  << block.path_ranges.size()
                  << std::endl;
        std::string consensus_id = (args::get(add_consensus) ? "Consensus." + std::to_string(block_id) : "");
        auto block_graph = smoothxg::smooth(graph, block, consensus_id);
        block_graph.to_gfa(std::cout);
    }
    */

    /*
    if (!args::get(xg_out).empty()) {
        std::ofstream out(args::get(xg_out));
        graph.serialize(out);
    }

    if (args::get(gfa_out)) {
        graph.to_gfa(std::cout);
    }
    */

    return 0;
}
