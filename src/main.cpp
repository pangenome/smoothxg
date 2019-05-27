/** \file version_main.cpp
 *
 * Defines the "vg version" subcommand, which evaluates graphs and alignments.
 */


#include <omp.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "args.hxx"
#include "sdsl/bit_vectors.hpp"
#include "xg.hpp"

using namespace std;
using namespace xg;

int main(int argc, char** argv) {
    args::ArgumentParser parser("xg: succinct static variation graph");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> gfa_in(parser, "FILE", "index the graph in this GFA file", {'g', "gfa-in"});
    args::ValueFlag<std::string> xg_out(parser, "FILE", "write the resulting xg index to this file", {'o', "out"});
    args::ValueFlag<std::string> xg_in(parser, "FILE", "read the xg index from this file", {'i', "in"});
    args::ValueFlag<std::string> base(parser, "BASE", "use this basename for temporary files during build", {'b', "base"});
    args::Flag gfa_out(parser, "FILE", "write the graph in GFA to stdout", {'G', "gfa-out"});
    args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
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

    //unique_ptr<XG> graph;
    XG graph;
    if (!args::get(xg_in).empty()) {
        std::ifstream in(args::get(xg_in));
        graph.load(in);
    } else if (!args::get(gfa_in).empty()) {
        graph.from_gfa(args::get(gfa_in));
    }

    if (!args::get(xg_out).empty()) {
        std::ofstream out(args::get(xg_out));
        graph.serialize(out);
    }

    if (args::get(gfa_out)) {
        graph.to_gfa(std::cout);
    }
    
    /*
    unique_ptr<XG> graph;
    //string file_name = argv[optind];
    if (in_name.empty()) assert(!vg_in.empty());
    if (vg_in == "-") {
        // Read VG from stdin
        graph = unique_ptr<XG>(new XG());
        graph->from_stream(std::cin, validate_graph, print_graph, store_threads, is_sorted_dag);
    } else if (vg_in.size()) {
        // Read VG from a file
        ifstream in;
        in.open(vg_in.c_str());
        graph = unique_ptr<XG>(new XG());
        graph->from_stream(in, validate_graph, print_graph, store_threads, is_sorted_dag);
    }

    if (in_name.size()) {
        get_input_file(in_name, [&](istream& in) {
            // Load from an XG file or - (stdin)
            graph = vg::io::VPKG::load_one<XG>(in);
        });
    }

    // Prepare structure tree for serialization
    unique_ptr<sdsl::structure_tree_node> structure;
    
    if (!report_name.empty()) {
        // We need to make a report, so we need the structure. Make a real tree
        // node. The unique_ptr handles deleting.
        structure = unique_ptr<sdsl::structure_tree_node>(new sdsl::structure_tree_node("name", "type"));
    }

    if(!vg_out.empty()) {
        if (graph.get() == nullptr) {
             cerr << "error [vg xg] no xg graph exists to convert; Try: vg xg -i graph.xg -X graph.vg" << endl;
             return 1;
        }
        
        VG converted;
        // Convert the xg graph to vg format
        convert_path_handle_graph(graph.get(), &converted);
        
        if (vg_out == "-") {
            converted.serialize_to_ostream(std::cout);
        } else {
            converted.serialize_to_file(vg_out);
        }
    }

    if (!out_name.empty()) {
        // Open a destination file if it is a file we want to write to
        ofstream out_file;
        if (out_name != "-") {
            out_file.open(out_name);
        }
        // Work out where to save to
        ostream& out = (out_name == "-") ? std::cout : out_file;
        
        // Encapsulate output in VPKG
        vg::io::VPKG::with_save_stream(out, "XG", [&](ostream& tagged) {
            // Serialize to the file while recording space usage to the structure.
            graph->serialize_and_measure(tagged, structure.get(), "xg");
        });
        
        out.flush();
    }

    if (!report_name.empty()) {
        // Save the report
        ofstream out;
        out.open(report_name.c_str());
        sdsl::write_structure_tree<HTML_FORMAT>(structure.get(), out, 0);
    }

    // queries
    if (node_sequence) {
        cout << node_id << ": " << graph->node_sequence(node_id) << endl;
    }
    if (!pos_for_char.empty()) {
        // extract the position from the string
        int64_t id;
        bool is_rev;
        size_t off;
        extract_pos(pos_for_char, id, is_rev, off);
        // then pick it up from the graph
        cout << graph->pos_char(id, is_rev, off) << endl;
    }
    if (!pos_for_substr.empty()) {
        int64_t id;
        bool is_rev;
        size_t off;
        size_t len;
        extract_pos_substr(pos_for_substr, id, is_rev, off, len);
        cout << graph->pos_substr(id, is_rev, off, len) << endl;
    }
    
    if (edges_from) {
        vector<Edge> edges = graph->edges_from(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << (edge.from_start()?"-":"+")
                 << " -> " << edge.to() << (edge.to_end()?"-":"+") << endl;
        }
    }
    if (edges_to) {
        vector<Edge> edges = graph->edges_to(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << (edge.from_start()?"-":"+")
                 << " -> " << edge.to() << (edge.to_end()?"-":"+") << endl;
        }
    }
    if (edges_of) {
        vector<Edge> edges = graph->edges_of(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << (edge.from_start()?"-":"+")
                 << " -> " << edge.to() << (edge.to_end()?"-":"+") << endl;
        }
    }
    if (edges_on_start) {
        vector<Edge> edges = graph->edges_on_start(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << (edge.from_start()?"-":"+")
                 << " -> " << edge.to() << (edge.to_end()?"-":"+") << endl;
        }
    }
    if (edges_on_end) {
        vector<Edge> edges = graph->edges_on_end(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << (edge.from_start()?"-":"+")
                 << " -> " << edge.to() << (edge.to_end()?"-":"+") << endl;
        }
    }

    if (node_context) {
        Graph g;
        graph->neighborhood(node_id, context_steps, g);
        if (text_output) {
            to_text(cout, g);
        } else {
            vector<Graph> gb = { g };
            vg::io::write_buffered(cout, gb, 0);
        }
    }

    if (!target.empty()) {
        string name;
        int64_t start, end;
        Graph g;
        parse_region(target, name, start, end);
        graph->get_path_range(name, start, end, g);
        graph->expand_context(g, context_steps);
        if (text_output) {
            to_text(cout, g);
        } else {
            vector<Graph> gb = { g };
            vg::io::write_buffered(cout, gb, 0);
        }
    }
    
    if (extract_threads) {
        list<XG::thread_t> threads;
        for (auto& p : graph->extract_threads(false)) {
            for (auto& t : p.second) {
                threads.push_back(t);
            }
        }
        for (auto& p : graph->extract_threads(true)) {
            for (auto& t : p.second) {
                threads.push_back(t);
            }
        }

        size_t thread_number = 0;
        for(XG::thread_t& thread : threads) {
            // Convert to a Path
            Path path;
            for(XG::ThreadMapping& m : thread) {
                // Convert all the mappings
                Mapping mapping;
                mapping.mutable_position()->set_node_id(m.node_id);
                mapping.mutable_position()->set_is_reverse(m.is_reverse);
                
                *(path.add_mapping()) = mapping;
            }
        
        
            // Give each thread a name
            path.set_name("_thread_" + to_string(thread_number++));
            
            // We need a Graph for serialization purposes. We do one chunk per
            // thread in case the threads are long.
            Graph g;
            
            *(g.add_path()) = path;
            
            // Dump the graph with its mappings. TODO: can we restrict these to
            // mappings to nodes we have already pulled out? Or pull out the
            // whole compressed graph?
            if (text_output) {
                to_text(cout, g);
            } else {
                vector<Graph> gb = { g };
                vg::io::write_buffered(cout, gb, 0);
            }
            
        }
    }

    if (!b_array_name.empty()) {
        // Dump B array
        ofstream out;
        out.open(b_array_name.c_str());
        graph->bs_dump(out);
    }
    */

    return 0;
}
