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
    args::ValueFlag<uint64_t> max_gap_length(parser, "N", "maximum gap length in chaining (default 1000)", {'l', "max-gap-length"});
    args::ValueFlag<double> max_mismatch_rate(parser, "FLOAT", "maximum allowed mismatch rate (default 0.1)", {'r', "max-mismatch-rate"});
    args::ValueFlag<double> chain_overlap(parser, "FLOAT", "maximum allowed query overlap between chains superchains (default 1.0 / disabled)", {'c', "chain-overlap-max"});
    args::ValueFlag<uint64_t> chain_min_anchors(parser, "N", "minimum number of anchors in a chain (3)", {'a', "chain-min-n-anchors"});
    args::ValueFlag<uint64_t> align_best_n(parser, "N", "align the best N superchains", {'n', "align-best-n"});
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

    uint64_t max_gap = args::get(max_gap_length)
        ? args::get(max_gap_length) : 1000;

    double mismatch_rate = args::get(max_mismatch_rate)
        ? args::get(max_mismatch_rate)
        : 0.2;

    double chain_overlap_max = args::get(chain_overlap)
        ? args::get(chain_overlap)
        : 1.0;

    uint64_t chain_min_n_anchors = args::get(chain_min_anchors)
        ? args::get(chain_min_anchors)
        : 3;

    uint64_t best_n = args::get(align_best_n) ? args::get(align_best_n) : 1;

    auto blocks = smoothxg::collinear_blocks(graph,
                                             max_gap,
                                             mismatch_rate,
                                             chain_min_n_anchors,
                                             chain_overlap_max);

    if (!args::get(xg_out).empty()) {
        std::ofstream out(args::get(xg_out));
        graph.serialize(out);
    }

    if (args::get(gfa_out)) {
        graph.to_gfa(std::cout);
    }

    return 0;
}
