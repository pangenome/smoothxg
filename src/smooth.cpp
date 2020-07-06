#include "smooth.hpp"


namespace smoothxg {


void smooth(const xg::XG& graph,
            const block_t& block) {

    // TODO we should take these as input
    std::int8_t poa_m = 5;
    std::int8_t poa_n = -4;
    std::int8_t poa_g = -8;
    std::int8_t poa_e = -6;
    std::int8_t poa_q = -10;
    std::int8_t poa_c = -4;
    std::uint8_t poa_algorithm = 0;
    
    auto poa_graph = spoa::createGraph();
    // collect sequences
    std::vector<std::string> seqs;
    std::vector<std::string> names;
    for (auto& path_range : block.path_ranges) {
        seqs.emplace_back();
        auto& seq = seqs.back();
        for (step_handle_t step = path_range.begin;
             step != path_range.end;
             step = graph.get_next_step(step)) {
            seq.append(graph.get_sequence(graph.get_handle_of_step(step)));
        }
        std::stringstream namess;
        namess << graph.get_path_name(graph.get_path_handle_of_step(path_range.begin))
               << "_" << graph.get_position_of_step(path_range.begin);
        names.push_back(namess.str());
    }
    // set up POA
    std::unique_ptr<spoa::AlignmentEngine> alignment_engine;
    try {
        alignment_engine = spoa::createAlignmentEngine(
            static_cast<spoa::AlignmentType>(poa_algorithm),
            poa_m, poa_n, poa_g, poa_e, poa_q, poa_c);
    } catch(std::invalid_argument& exception) {
        std::cerr << exception.what() << std::endl;
        assert(false);
    }
    // run POA
    std::size_t max_sequence_size = 0;
    for (auto& seq : seqs) {
        max_sequence_size = std::max(max_sequence_size, seq.size());
    }
    alignment_engine->prealloc(max_sequence_size, 4);
    for (auto& seq : seqs) {
        auto alignment = alignment_engine->align(seq, poa_graph);
        try {
            poa_graph->add_alignment(alignment, seq); // could give weight
        } catch(std::invalid_argument& exception) {
            std::cerr << exception.what() << std::endl;
            assert(false);
        }
    }
    // force consensus genertion for graph annotation
    std::string consensus = poa_graph->generate_consensus();
    // write the graph, with consensus as a path if requested
    poa_graph->print_gfa(std::cout, names, true);
}

}
