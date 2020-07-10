#include "smooth.hpp"


namespace smoothxg {

void smooth_and_lace(const xg::XG& graph,
                     const std::vector<block_t>& blocks) {
    // for each block, smooth into a file
    // record the start and end points of all the path ranges and the consensus
    // 
}

void smooth(const xg::XG& graph,
            const block_t& block,
            std::ostream& out) {

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
        // TODO determine alignment orientation somehow!!
        auto alignment = alignment_engine->align(seq, poa_graph);
        try {
            poa_graph->add_alignment(alignment, seq); // could give weight
        } catch(std::invalid_argument& exception) {
            std::cerr << exception.what() << std::endl;
            assert(false);
        }
    }
    // todo make the consensus generation optional
    
    // force consensus genertion for graph annotation
    std::string consensus = poa_graph->generate_consensus();
    // write the graph, with consensus as a path if requested
    //poa_graph->print_gfa(poa_graph, std::cout, names, true);
    // optionally write in a different format?
    // or build a graph?
    //write_gfa(poa_graph, out, names, true);
    odgi::graph_t output_graph;
    build_odgi(poa_graph, output_graph, names, true);
    odgi::algorithms::unchop(output_graph);
    output_graph.apply_ordering(odgi::algorithms::topological_order(&output_graph), true);
    output_graph.to_gfa(out);
}

void write_gfa(std::unique_ptr<spoa::Graph>& graph,
               std::ostream& out,
               const std::vector<std::string>& sequence_names,
               bool include_consensus) {

    auto& nodes = graph->nodes();
    std::vector<std::int32_t> in_consensus(nodes.size(), -1);
    std::int32_t rank = 0;
    auto consensus = graph->consensus();
    for (const auto& id: consensus) {
        in_consensus[id] = rank++;
    }

    out << "H" << "\t" << "VN:Z:1.0" << std::endl;

    for (std::uint32_t i = 0; i < nodes.size(); ++i) {
        out << "S" << "\t" << i+1 << "\t" << static_cast<char>(graph->decoder(nodes[i]->code()));
        if (in_consensus[i] != -1) {
            out << "\t" << "ic:Z:true";
        }
        out << std::endl;
        for (const auto& edge: nodes[i]->out_edges()) {
            out << "L" << "\t" << i+1 << "\t" << "+" << "\t" << edge->end_node_id()+1 << "\t" << "+" << "\t" << "0M" << "\t"
                << "ew:f:" << edge->total_weight();
            if (in_consensus[i] + 1 == in_consensus[edge->end_node_id()]) {
                out << "\t" << "ic:Z:true";
            }
            out << std::endl;
        }
    }

    for (std::uint32_t i = 0; i < sequence_names.size(); ++i) {
        out << "P" << "\t" << sequence_names[i] << "\t";
        std::uint32_t node_id = graph->sequences_begin_nodes_ids()[i];
        while (true) {
            out << node_id+1 << "+";
            if (!nodes[node_id]->successor(node_id, i)) {
                break;
            } else {
                out << ",";
            }
        }
        out << "\t" << "*" << std::endl;
    }

    if (include_consensus) {
        out << "P" << "\t" << "Consensus" << "\t";
        for (auto& id : graph->consensus()) {
            out << id+1 << "+";
        }
        out << "\t" << "*" << std::endl;
    }
}

void build_odgi(std::unique_ptr<spoa::Graph>& graph,
                odgi::graph_t& output,
                const std::vector<std::string>& sequence_names,
                bool include_consensus) {

    auto& nodes = graph->nodes();
    std::vector<std::int32_t> in_consensus(nodes.size(), -1);
    std::int32_t rank = 0;
    auto consensus = graph->consensus();
    for (const auto& id: consensus) {
        in_consensus[id] = rank++;
    }

    for (std::uint32_t i = 0; i < nodes.size(); ++i) {
        std::string seq = std::string(1, static_cast<char>(graph->decoder(nodes[i]->code())));
        output.create_handle(seq, i+1);
    }

    for (std::uint32_t i = 0; i < nodes.size(); ++i) {
        for (const auto& edge: nodes[i]->out_edges()) {
            output.create_edge(output.get_handle(i+1), output.get_handle(edge->end_node_id()+1));
        }
    }

    for (std::uint32_t i = 0; i < sequence_names.size(); ++i) {
        path_handle_t p = output.create_path_handle(sequence_names[i]);
        std::uint32_t node_id = graph->sequences_begin_nodes_ids()[i];
        while (true) {
            output.append_step(p, output.get_handle(node_id+1));
            if (!nodes[node_id]->successor(node_id, i)) {
                break;
            }
        }
    }

    if (include_consensus) {
        path_handle_t p = output.create_path_handle("Consensus"); // TODO configure name
        for (auto& id : graph->consensus()) {
            output.append_step(p, output.get_handle(id+1));
        }
    }
}


}
