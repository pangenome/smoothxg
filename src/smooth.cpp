#include "smooth.hpp"


namespace smoothxg {

odgi::graph_t smooth(const xg::XG& graph,
                     const block_t& block,
                     const std::string& consensus_name) {

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
        // TODO
        // preserve mapping between these entities and the global path handles
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
    //std::vector<bool> aln_is_reverse;
    for (auto& seq : seqs) {
        // TODO determine alignment orientation somehow!!!!!!!!
        // or try both orientations here
        // we'll need to record the orientation in the path somehow
        // to structure the lacing
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
    // write the graph, with consensus as a path
    odgi::graph_t output_graph;
    // convert the poa graph into our output format
    build_odgi(poa_graph, output_graph, names, consensus_name, !consensus_name.empty());
    // normalize the representation, allowing for nodes > 1bp
    odgi::algorithms::unchop(output_graph);
    // order the graph
    output_graph.apply_ordering(odgi::algorithms::topological_order(&output_graph), true);
    //output_graph.to_gfa(out);
    return output_graph;
}

odgi::graph_t smooth_and_lace(const xg::XG& graph,
                              const std::vector<block_t>& blocks,
                              const std::string& consensus_base_name) {
    //
    // record the start and end points of all the path ranges and the consensus
    //
    std::vector<odgi::graph_t> block_graphs;
    block_graphs.reserve(blocks.size());
    std::vector<path_position_range_t> path_mapping;
    std::vector<path_position_range_t> consensus_mapping;
    bool add_consensus = !consensus_base_name.empty();
    uint64_t block_id = 0;
    for (auto& block : blocks) {
        std::string consensus_name = consensus_base_name + std::to_string(block_id);
        block_graphs.push_back(smooth(graph, block, consensus_name));
        auto& block_graph = block_graphs.back();
        // record the start and end paths
        // nb: the path order is the same in the input block and output graph
        uint64_t path_id = 0;
        for (auto& path_range : block.path_ranges) {
            auto path_handle = graph.get_path_handle_of_step(path_range.begin);
            auto last_step = graph.get_previous_step(path_range.end);
            path_mapping.push_back({
                    path_handle, // target path
                    graph.get_position_of_step(path_range.begin), // start position
                    (graph.get_position_of_step(last_step) // end position
                     + graph.get_length(graph.get_handle_of_step(last_step))),
                    as_path_handle(++path_id),
                    block_id
                });
        }
        // make the graph
        
        // record the consensus path
        if (add_consensus) {
            auto consensus_handle = block_graph.get_path_handle(consensus_name);
            uint64_t path_end = 0;
            block_graph.for_each_step_in_path(
                consensus_handle,
                [&](const step_handle_t& step) {
                    path_end += block_graph.get_length(block_graph.get_handle_of_step(step));
                });
            consensus_mapping.push_back({
                    as_path_handle(0), // consensus = 0 path handle
                    0, // start position
                    path_end, // end position
                    consensus_handle,
                    block_id
                });
        }
        // increment our block id
        ++block_id;
    }
    // sort the path range mappings by path handle id, then start position
    // this will allow us to walk through them in order
    ips4o::parallel::sort(
        path_mapping.begin(),
        path_mapping.end(),
        [](const path_position_range_t& a,
           const path_position_range_t& b) {
            auto& a_id = as_integer(a.base_path);
            auto& b_id = as_integer(b.base_path);
            return (a_id < b_id || a_id == b_id && a.start_pos < b.start_pos);
        });
    // build the sequence and edges into the output graph
    odgi::graph_t smoothed;
    // add the nodes and edges to the graph
    std::vector<uint64_t> id_mapping;
    for (auto& block : block_graphs) {
        // record the id translation
        uint64_t id_trans = smoothed.get_node_count();
        id_mapping.push_back(id_trans);
        //std::cerr << "block graph " << block.get_node_count() << std::endl;
        block.for_each_handle(
            [&](const handle_t& h) {
                smoothed.create_handle(block.get_sequence(h));
            });
        block.for_each_edge(
            [&](const edge_t& e) {
                smoothed.create_edge(
                    smoothed.get_handle(id_trans + block.get_id(e.first)),
                    smoothed.get_handle(id_trans + block.get_id(e.second))
                    );
            });
    }
    // then for each path, ensure that it's embedded in the graph by walking through its block segments in order
    // and linking them up in the output graph
    for (uint64_t i = 0; i < path_mapping.size(); ++i) {
        path_position_range_t* pos_range = &path_mapping[i];
        path_position_range_t* last_pos_range = nullptr;
        // add the path to the graph
        // walk the path from start to end
        do {
            // if we find a segment that's not included in any path, we'll add it to the final graph and link it in
            if (last_pos_range != nullptr) {
                // if we have a gap in length, collect the sequence in the gap and add it to the graph as a node
                // then add it as a traversal to the path
            }
            // write the path steps into the graph using the id translation
            
            last_pos_range = pos_range;
        } while (path_mapping[++i].base_path == pos_range->base_path);
    }
    return smoothed;
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
                const std::string& consensus_name,
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
        path_handle_t p = output.create_path_handle(consensus_name);
        for (auto& id : graph->consensus()) {
            output.append_step(p, output.get_handle(id+1));
        }
    }
}


}
