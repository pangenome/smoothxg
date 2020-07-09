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

using namespace std;
using namespace xg;

int main(int argc, char** argv) {
    args::ArgumentParser parser("smoothxg: collinear block finder and graph consensus generator");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> gfa_in(parser, "FILE", "index the graph in this GFA file", {'g', "gfa-in"});
    args::ValueFlag<std::string> xg_out(parser, "FILE", "write the resulting xg index to this file", {'o', "out"});
    args::ValueFlag<std::string> xg_in(parser, "FILE", "read the xg index from this file", {'i', "in"});
    args::ValueFlag<std::string> base(parser, "BASE", "use this basename for temporary files during build", {'b', "base"});
    args::Flag gfa_out(parser, "FILE", "write the graph in GFA to stdout", {'G', "gfa-out"});
    args::ValueFlag<uint64_t> _max_block_weight(parser, "N", "maximum seed sequence in block (default: 10000)", {'w', "max-block-weight"});
    args::ValueFlag<uint64_t> _max_block_jump(parser, "N", "maximum path jump to include in block (default: 1000)", {'j', "max-path-jump"});
    args::ValueFlag<uint64_t> _min_subpath(parser, "N", "minimum length of a subpath to include in partial order alignment (default: 16)", {'k', "min-subpath"});
    args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    args::Flag validate(parser, "validate", "validate construction", {'V', "validate"});
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
        graph.from_gfa(args::get(gfa_in), args::get(validate),
                       args::get(base).empty() ? args::get(gfa_in) : args::get(base));
    }

    uint64_t max_block_weight = args::get(_max_block_weight) ? args::get(_max_block_weight) : 10000;
    uint64_t max_block_jump = args::get(_max_block_jump) ? args::get(_max_block_jump) : 1000;
    uint64_t min_subpath = args::get(_min_subpath) ? args::get(_min_subpath) : 16;

    auto blocks = smoothxg::smoothable_blocks(graph,
                                              max_block_weight,
                                              max_block_jump);

    uint64_t block_id = 0;
    for (auto& block : blocks) {
        std::cout << "block" << block_id++ << "\t"
                  << block.total_path_length << "\t"
                  << block.max_path_length << "\t"
                  << graph.get_id(block.handles.front())
                  << "-" << graph.get_id(block.handles.back()) << "\t"
                  << block.path_ranges.size()
                  << std::endl;
        smoothxg::smooth(graph, block, std::cout);
    }

    if (!args::get(xg_out).empty()) {
        std::ofstream out(args::get(xg_out));
        graph.serialize(out);
    }

    if (args::get(gfa_out)) {
        graph.to_gfa(std::cout);
    }

    return 0;
}
