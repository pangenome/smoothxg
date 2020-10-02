#include "smooth.hpp"
#include <cstring>

#include "deps/abPOA/src/seq.h"

namespace smoothxg {

// klib stuff copied from abpoa_graph.c
KDQ_INIT(int)
#define kdq_int_t kdq_t(int)

// to write each block to a FASTA and TSV
//#define SMOOTH_WRITE_BLOCKS_FASTA true

static inline int ilog2_64(abpoa_para_t *abpt, uint64_t v) {
    uint64_t t, tt;
    if ((tt = v >> 32)) {
        return (t = tt >> 16) ? 48 + abpt->LogTable65536[t]
                              : 32 + abpt->LogTable65536[tt];
    }
    return (t = v >> 16) ? 16 + abpt->LogTable65536[t] : abpt->LogTable65536[v];
}

odgi::graph_t smooth_abpoa(const xg::XG &graph, const block_t &block, const uint64_t &block_id,
                           int poa_m, int poa_n, int poa_g,
                           int poa_e, int poa_q, int poa_c,
                           bool local_alignment,
                           const std::string &consensus_name) {

    // collect sequences
    std::vector<std::string> seqs;
    std::vector<std::string> names;
    for (auto &path_range : block.path_ranges) {
        seqs.emplace_back();
        auto &seq = seqs.back();
        for (step_handle_t step = path_range.begin; step != path_range.end;
             step = graph.get_next_step(step)) {
            seq.append(graph.get_sequence(graph.get_handle_of_step(step)));
        }
        std::stringstream namess;
        namess << graph.get_path_name(
                      graph.get_path_handle_of_step(path_range.begin))
               << "_" << graph.get_position_of_step(path_range.begin);
        names.push_back(namess.str());
    }

#ifdef SMOOTH_WRITE_BLOCKS_FASTA
    {
        std::string s = "smoothxg_block_" + std::to_string(block_id) + ".fa";
        std::ofstream fasta(s.c_str());
        for (uint64_t i = 0; i < seqs.size(); ++i) {
            fasta << ">" << names[i] << " " << seqs[i].size() << std::endl
                  << seqs[i] << std::endl;
        }
        fasta.close();
        const VectorizableHandleGraph &vec_graph =
            dynamic_cast<const VectorizableHandleGraph &>(graph);
        std::string v = "smoothxg_block_" + std::to_string(block_id) + ".tsv";
        std::ofstream vs(v.c_str());
        vs << "path.name\tstep.rank\tpos\tnode.id\tnode.pos\tvisit"
           << std::endl;
        for (auto &path_range : block.path_ranges) {
            std::string path_name = graph.get_path_name(
                graph.get_path_handle_of_step(path_range.begin));
            uint64_t rank = 0;
            uint64_t pos = 0;
            std::map<uint64_t, uint64_t> visits;
            for (step_handle_t step = path_range.begin; step != path_range.end;
                 step = graph.get_next_step(step)) {
                handle_t h = graph.get_handle_of_step(step);
                uint64_t id = graph.get_id(h);
                int64_t node_pos = vec_graph.node_vector_offset(id);
                auto &visit = visits[id];
                vs << path_name << "\t" << rank++ << "\t" << pos << "\t" << id
                   << "\t" << node_pos << "\t" << visit << std::endl;
                ++visit;
                pos += graph.get_length(graph.get_handle_of_step(step));
            }
        }
        vs.close();
    }
#endif

    // set up POA
    // done...
    // run POA
    std::size_t max_sequence_size = 0;
    for (auto &seq : seqs) {
        max_sequence_size = std::max(max_sequence_size, seq.size());
    }

    odgi::graph_t output_graph;
    // if the graph would be empty, bail out
    if (max_sequence_size == 0) {
        return output_graph;
    }

    bool generate_consensus = !consensus_name.empty();

    // initialize abPOA
    abpoa_t *ab = abpoa_init();
    // initialize abPOA parameters
    abpoa_para_t *abpt = abpoa_init_para();
    // if we want to do local alignments
    if (local_alignment) abpt->align_mode = ABPOA_LOCAL_MODE;
    //abpt->zdrop = 100; // could be useful in local mode
    //abpt->end_bonus = 100; // also useful in local mode
    abpt->rev_cigar = 0;
    abpt->out_gfa = 1; // must be set to get the graph
    //abpt->out_msa = 1; // must be set when we extract the MSA
    abpt->out_cons = generate_consensus;
    abpt->amb_strand = 1;
    abpt->match = poa_m;
    abpt->mismatch = poa_n;
    abpt->gap_open1 = poa_g;
    abpt->gap_open2 = poa_q;
    abpt->gap_ext1 = poa_e;
    abpt->gap_ext2 = poa_c;

    // finalize parameters
    abpoa_post_set_para(abpt);

    std::vector<char *> seqs_;
    // transform so that we have an interface between C++ and C
    std::transform(seqs.begin(), seqs.end(), std::back_inserter(seqs_),
                   [](const std::string& s) { return (char*)s.c_str();});

    // collect sequence length, transform ACGT to 0123
    int n_seqs = seqs.size();
    int *seq_lens = (int *)malloc(sizeof(int) * n_seqs);
    uint8_t **bseqs = (uint8_t **)malloc(sizeof(uint8_t *) * n_seqs);
    for (int i = 0; i < n_seqs; ++i) {
        seq_lens[i] = seqs[i].size();
        bseqs[i] = (uint8_t *)malloc(sizeof(uint8_t) * seq_lens[i]);
        for (int j = 0; j < seq_lens[i]; ++j) {
            bseqs[i][j] =
                nst_nt4_table[(int)seqs_[i][j]]; // TODO we make a c_str for
                                                 // every char in the string
        }
    }
    uint8_t **cons_seq;
    int **cons_cov, *cons_l, cons_n = 0;
    uint8_t **msa_seq;
    int msa_l = 0;
    // perform abpoa-msa
    int i, tot_n = n_seqs;
    uint8_t *is_rc = (uint8_t *)_err_malloc(n_seqs * sizeof(uint8_t));
    abpoa_reset_graph(ab, abpt, seq_lens[0]);
    
    std::vector<bool> aln_is_reverse;
    for (i = 0; i < n_seqs; ++i) {
        abpoa_res_t res;
        res.graph_cigar = 0, res.n_cigar = 0, res.is_rc = 0;
        abpt->rev_cigar = 0;
        abpoa_align_sequence_to_graph(ab, abpt, bseqs[i], seq_lens[i], &res);
        abpoa_add_graph_alignment(ab, abpt, bseqs[i], seq_lens[i], res, i,
                                  n_seqs);
        is_rc[i] = res.is_rc;
        if (res.is_rc) {
            aln_is_reverse.push_back(true);
            //std::cerr << "is_rc" << std::endl;
        } else {
            aln_is_reverse.push_back(false);
            // std::cerr << "is_rc_not" << std::endl;
        }
        if (res.n_cigar) {
            free(res.graph_cigar);
        }
    }

    if (generate_consensus) {
        abpoa_generate_consensus(ab, abpt, tot_n, NULL, NULL, NULL, NULL, NULL);
        if (ab->abg->is_called_cons == 0) {
            // TODO Abort mission here?
            err_printf("ERROR: no consensus sequence generated.\n");
            exit(1);
        }
        aln_is_reverse.push_back(false);
    }
    free(is_rc);

    /*
    fprintf(stdout, "=== output to variables ===\n");
    for (int i = 0; i < cons_n; ++i) {
        fprintf(stdout, ">Consensus_sequence\n");
        for (int j = 0; j < cons_l[i]; ++j)
            fprintf(stdout, "%c", nst_nt256_table[cons_seq[i][j]]);
        fprintf(stdout, "\n");
    }

    fprintf(stdout, ">Multiple_sequence_alignment\n");
    for (i = 0; i < n_seqs; ++i) {
        for (int j = 0; j < msa_l; ++j) {
            fprintf(stdout, "%c", "ACGTN-"[msa_seq[i][j]]);
        }
        fprintf(stdout, "\n");
    }
    */

    if (cons_n) {
        for (i = 0; i < cons_n; ++i) {
            free(cons_seq[i]);
            free(cons_cov[i]);
        }
        free(cons_seq);
        free(cons_cov);
        free(cons_l);
    }
    if (msa_l) {
        for (i = 0; i < n_seqs; ++i) {
            free(msa_seq[i]);
        }
        free(msa_seq);
    }

    for (i = 0; i < n_seqs; ++i) {
        free(bseqs[i]);
    }
    free(bseqs);
    free(seq_lens);

    build_odgi_abPOA(ab, abpt, output_graph, names, aln_is_reverse,
                     consensus_name, generate_consensus);

    abpoa_free(ab, abpt);
    abpoa_free_para(abpt);
    // normalize the representation, allowing for nodes > 1bp
    odgi::algorithms::unchop(output_graph);
    // order the graph
    output_graph.apply_ordering(
        odgi::algorithms::topological_order(&output_graph), true);

    // output_graph.to_gfa(std::cout);
    return output_graph;
}

odgi::graph_t smooth_spoa(const xg::XG &graph, const block_t &block,
                          const uint64_t &block_id,
                          std::int8_t poa_m, std::int8_t poa_n, std::int8_t poa_g,
                          std::int8_t poa_e, std::int8_t poa_q, std::int8_t poa_c,
                          bool local_alignment,
                          const std::string &consensus_name) {

    std::uint8_t spoa_algorithm = local_alignment ? 0 : 1;
    std::unique_ptr<spoa::AlignmentEngine> alignment_engine
        = spoa::createAlignmentEngine(
            static_cast<spoa::AlignmentType>(spoa_algorithm),
            poa_m, poa_n, poa_g, poa_e, poa_q, poa_c);

    auto poa_graph = spoa::createGraph();
    // collect sequences
    std::vector<std::string> seqs;
    std::vector<std::string> names;
    for (auto &path_range : block.path_ranges) {
        seqs.emplace_back();
        auto &seq = seqs.back();
        for (step_handle_t step = path_range.begin; step != path_range.end;
             step = graph.get_next_step(step)) {
            seq.append(graph.get_sequence(graph.get_handle_of_step(step)));
        }
        std::stringstream namess;
        namess << graph.get_path_name(
                      graph.get_path_handle_of_step(path_range.begin))
               << "_" << graph.get_position_of_step(path_range.begin);
        names.push_back(namess.str());
    }
    /*
    std::string s = "smoothxg_block_" + std::to_string(block_id) + ".fa";
    std::ofstream fasta(s.c_str());
    for (uint64_t i = 0; i < seqs.size(); ++i) {
        fasta << ">" << names[i] << " " << seqs[i].size() << std::endl
              << seqs[i] << std::endl;
    }
    fasta.close();
    */
    // set up POA
    // done...
    // run POA
    std::size_t max_sequence_size = 0;
    for (auto &seq : seqs) {
        max_sequence_size = std::max(max_sequence_size, seq.size());
    }
    odgi::graph_t output_graph;
    // if the graph would be empty, bail out
    if (max_sequence_size == 0) {
        return output_graph;
    }
    // preallocation does not seem to help, and it consumes a lot of memory
    // relative to progressive allocation
    // alignment_engine->prealloc(max_sequence_size, 4);
    std::vector<bool> aln_is_reverse;
    int i = 0;
    for (auto &seq : seqs) {
        // std::cerr << names[i++] << "\t" << seq << std::endl;
        // TODO determine alignment orientation somehow!!!!!!!!
        // or try both orientations here
        // we'll need to record the orientation in the path somehow
        // to structure the lacing
        std::int32_t score_fwd = 0;
        auto alignment_fwd =
            alignment_engine->align(seq, poa_graph, &score_fwd);
        auto rev_seq = odgi::reverse_complement(seq);
        std::int32_t score_rev = 0;
        auto alignment_rev =
            alignment_engine->align(rev_seq, poa_graph, &score_rev);
        try {
            // could give weight here to influence consensus
            if (score_fwd >= score_rev) {
                poa_graph->add_alignment(alignment_fwd, seq);
                aln_is_reverse.push_back(false);
            } else {
                poa_graph->add_alignment(alignment_rev, rev_seq);
                aln_is_reverse.push_back(true);
            }
        } catch (std::invalid_argument &exception) {
            std::cerr << exception.what() << std::endl;
            assert(false);
        }
    }
    // todo make the consensus generation optional
    // ...
    // force consensus genertion for graph annotation
    std::string consensus = poa_graph->generate_consensus();
    aln_is_reverse.push_back(false);
    // write the graph, with consensus as a path
    // odgi::graph_t output_graph;
    // convert the poa graph into our output format
    // poa_graph->print_gfa(std::cout, names, true);
    build_odgi(poa_graph, output_graph, names, aln_is_reverse, consensus_name,
               !consensus_name.empty());
    // normalize the representation, allowing for nodes > 1bp
    odgi::algorithms::unchop(output_graph);
    // order the graph
    output_graph.apply_ordering(
        odgi::algorithms::topological_order(&output_graph), true);
    // output_graph.to_gfa(out);
    return output_graph;
}

odgi::graph_t smooth_and_lace(const xg::XG &graph,
                              const std::vector<block_t> &blocks,
                              int poa_m, int poa_n,
                              int poa_g, int poa_e,
                              int poa_q, int poa_c,
                              bool local_alignment,
                              bool use_abpoa,
                              const std::string &consensus_base_name) {

    //
    // record the start and end points of all the path ranges and the consensus
    //
    std::vector<odgi::graph_t> block_graphs;
    block_graphs.resize(blocks.size());
    std::vector<path_position_range_t> path_mapping;
    std::vector<path_position_range_t> consensus_mapping;
    bool add_consensus = !consensus_base_name.empty();
    std::mutex path_mapping_mutex, consensus_mapping_mutex, logging_mutex;
    uint64_t thread_count = odgi::get_thread_count();

    paryfor::parallel_for<uint64_t>(
        0, blocks.size(), thread_count, [&](uint64_t block_id, int tid) {
            auto &block = blocks[block_id];

            { // if (block_id % 100 == 0) {
                std::lock_guard<std::mutex> guard(logging_mutex);
                std::cerr
                    << "[smoothxg::smooth_and_lace] applying " << (use_abpoa ? "abPOA" : "SPOA")
                    << " (" << (local_alignment ? "local" : "global") << " alignment mode)"
                    << " to block " << block_id << "/" << blocks.size() << " " << std::fixed
                    << std::showpoint << std::setprecision(3)
                    << (float)block_id / (float)blocks.size() * 100 << "%\r";
            }

            std::string consensus_name =
                consensus_base_name + std::to_string(block_id);
            // std::cerr << "on block " << block_id+1 << " of " << blocks.size()
            // << std::endl;
            auto &block_graph = block_graphs[block_id];

            if (use_abpoa) {
                block_graph = smooth_abpoa(graph,
                                           block,
                                           block_id,
                                           poa_m,
                                           poa_n,
                                           poa_g,
                                           poa_e,
                                           poa_q,
                                           poa_c,
                                           local_alignment,
                                           consensus_name);
            } else {
                block_graph = smooth_spoa(graph,
                                          block,
                                          block_id,
                                          poa_m,
                                          -poa_n,
                                          -poa_g,
                                          -poa_e,
                                          -poa_q,
                                          -poa_c,
                                          local_alignment,
                                          consensus_name);
            }

            // std::cerr << std::endl;
            // std::cerr << "After block graph. Exiting for now....." <<
            // std::endl; exit(0);
            if (block_graph.get_node_count() > 0) {
                // auto& block_graph = block_graphs.back();
                // record the start and end paths
                // nb: the path order is the same in the input block and output
                // graph
                uint64_t path_id = 0;
                for (auto &path_range : block.path_ranges) {
                    auto path_handle =
                        graph.get_path_handle_of_step(path_range.begin);
                    auto last_step = graph.get_previous_step(path_range.end);
                    {
                        std::lock_guard<std::mutex> guard(path_mapping_mutex);
                        path_mapping.push_back(
                            {path_handle, // target path
                             graph.get_position_of_step(
                                 path_range.begin), // start position
                             (graph.get_position_of_step(
                                  last_step) // end position
                              + graph.get_length(
                                    graph.get_handle_of_step(last_step))),
                             path_range.begin, path_range.end,
                             as_path_handle(++path_id), block_id});
                    }
                }
                // make the graph

                // record the consensus path
                if (add_consensus) {
                    auto consensus_handle =
                        block_graph.get_path_handle(consensus_name);
                    uint64_t path_end = 0;
                    step_handle_t empty_step;
                    as_integers(empty_step)[0] = 0;
                    as_integers(empty_step)[1] = 0;
                    block_graph.for_each_step_in_path(
                        consensus_handle, [&](const step_handle_t &step) {
                            path_end += block_graph.get_length(
                                block_graph.get_handle_of_step(step));
                        });
                    {
                        std::lock_guard<std::mutex> guard(
                            consensus_mapping_mutex);
                        consensus_mapping.push_back(
                            {as_path_handle(0), // consensus = 0 path handle
                             0,                 // start position
                             path_end,          // end position
                             empty_step, empty_step, consensus_handle,
                             block_id});
                    }
                }
            }
        });

    std::cerr << "[smoothxg::smooth_and_lace] applying " << (use_abpoa ? "abPOA" : "SPOA")
              << " (" << (local_alignment ? "local" : "global") << " alignment mode)"
              << " to block " << blocks.size() << "/" << blocks.size() << " " << std::fixed
              << std::showpoint << std::setprecision(3) << 100.0 << "%"
              << std::endl;

    std::cerr << "[smoothxg::smooth_and_lace] sorting path_mappings"
              << std::endl;
    // sort the path range mappings by path handle id, then start position
    // this will allow us to walk through them in order
    ips4o::parallel::sort(
        path_mapping.begin(), path_mapping.end(),
        [](const path_position_range_t &a, const path_position_range_t &b) {
            auto &a_id = as_integer(a.base_path);
            auto &b_id = as_integer(b.base_path);
            return (a_id < b_id || a_id == b_id && a.start_pos < b.start_pos);
        });
    // build the sequence and edges into the output graph
    odgi::graph_t smoothed;
    // add the nodes and edges to the graph
    std::vector<uint64_t> id_mapping;
    std::cerr << "[smoothxg::smooth_and_lace] building final graph"
              << std::endl;
    uint64_t j = 0;
    for (auto &block : block_graphs) {
        uint64_t id_trans = smoothed.get_node_count();
        { // if (j % 100 == 0) {
            std::cerr << "[smoothxg::smooth_and_lace] adding graph " << j << "/"
                      << block_graphs.size() << " " << std::fixed
                      << std::showpoint << std::setprecision(3)
                      << (float)j / (float)block_graphs.size() << "%\r";
        }
        ++j;
        // record the id translation
        id_mapping.push_back(id_trans);
        if (block.get_node_count() == 0) {
            continue;
        }
        block.for_each_handle([&](const handle_t &h) {
            smoothed.create_handle(block.get_sequence(h));
        });
        block.for_each_edge([&](const edge_t &e) {
            smoothed.create_edge(
                smoothed.get_handle(id_trans + block.get_id(e.first)),
                smoothed.get_handle(id_trans + block.get_id(e.second)));
        });
    }
    std::cerr << "[smoothxg::smooth_and_lace] adding graph " << j++ << "/"
              << block_graphs.size() << " 100.000%" << std::endl;
    // then for each path, ensure that it's embedded in the graph by walking
    // through its block segments in order and linking them up in the output
    // graph
    for (uint64_t i = 0; i < path_mapping.size(); ++i) {
        {
            std::cerr << "[smoothxg::smooth_and_lace] embedding path fragment "
                      << i << "/" << path_mapping.size() << "\r";
        }
        path_position_range_t *pos_range = &path_mapping[i];
        path_position_range_t *last_pos_range = nullptr;
        step_handle_t last_step = {0, 0};
        uint64_t last_end_pos = 0;
        // add the path to the graph
        path_handle_t smoothed_path = smoothed.create_path_handle(
            graph.get_path_name(pos_range->base_path));
        // walk the path from start to end
        while (true) {
            // if we find a segment that's not included in any block, we'll add
            // it to the final graph and link it in to do so, we detect a gap in
            // length, collect the sequence in the gap and add it to the graph
            // as a node then add it as a traversal to the path
            if (pos_range->start_pos - last_end_pos > 0) {
                step_handle_t last;
                if (last_pos_range != nullptr) {
                    // we can iterate from the last path end
                    last = last_pos_range->end_step;
                } else {
                    // iterate from the path start to here
                    last = graph.path_begin(pos_range->base_path);
                }
                // 1) collect sequence
                std::string seq;
                for (step_handle_t step = last; step != pos_range->start_step;
                     step = graph.get_next_step(step)) {
                    seq.append(
                        graph.get_sequence(graph.get_handle_of_step(step)));
                }
                // 2) create node
                handle_t h = smoothed.create_handle(seq);
                // 3) append to path in smoothed
                smoothed.append_step(smoothed_path, h);
                if (as_integers(last_step)[0] != 0) {
                    smoothed.create_edge(smoothed.get_handle_of_step(last_step),
                                         h);
                }
                last_step = smoothed.path_back(smoothed_path);
            }
            // write the path steps into the graph using the id translation
            auto &block = block_graphs[pos_range->target_graph_id];
            auto &id_trans = id_mapping[pos_range->target_graph_id];
            bool first = true;
            block.for_each_step_in_path(
                pos_range->target_path, [&](const step_handle_t &step) {
                    handle_t h = block.get_handle_of_step(step);
                    handle_t t = smoothed.get_handle(block.get_id(h) + id_trans,
                                                     block.get_is_reverse(h));
                    smoothed.append_step(smoothed_path, t);
                    if (first) {
                        first = false;
                        // create edge between last and curr
                        if (as_integers(last_step)[0] != 0) {
                            smoothed.create_edge(
                                smoothed.get_handle_of_step(last_step), t);
                        }
                    }
                });
            last_step = smoothed.path_back(smoothed_path);
            last_pos_range = pos_range;
            last_end_pos = pos_range->end_pos;
            if (i + 1 == path_mapping.size() ||
                path_mapping[i + 1].base_path != pos_range->base_path) {
                break;
            } else {
                ++i;
                pos_range = &path_mapping[i];
            }
        }
        // now add in any final sequence in the path
        // and add it to the path, add the edge
        if (graph.get_path_length(pos_range->base_path) > last_end_pos) {
            std::string seq;
            for (step_handle_t step = pos_range->end_step;
                 step != graph.path_end(pos_range->base_path);
                 step = graph.get_next_step(step)) {
                seq.append(graph.get_sequence(graph.get_handle_of_step(step)));
            }
            handle_t h = smoothed.create_handle(seq);
            smoothed.create_edge(smoothed.get_handle_of_step(last_step), h);
            smoothed.append_step(smoothed_path, h);
        }
    }
    std::cerr << "[smoothxg::smooth_and_lace] embedding path fragment "
              << path_mapping.size() << "/" << path_mapping.size() << std::endl;
    // now verify that smoothed has paths that are equal to the base graph
    // and that all the paths are fully embedded in the graph
    std::cerr << "[smoothxg::smooth_and_lace] verifying paths" << std::endl;
    smoothed.for_each_path_handle([&](const path_handle_t &path) {
        // collect sequence
        std::string orig_seq, smoothed_seq;
        graph.for_each_step_in_path(
            graph.get_path_handle(smoothed.get_path_name(path)),
            [&](const step_handle_t &step) {
                orig_seq.append(
                    graph.get_sequence(graph.get_handle_of_step(step)));
            });
        smoothed.for_each_step_in_path(path, [&](const step_handle_t &step) {
            smoothed_seq.append(
                smoothed.get_sequence(smoothed.get_handle_of_step(step)));
        });
        if (orig_seq != smoothed_seq) {
            std::cerr << "[smoothxg] error! path "
                      << smoothed.get_path_name(path)
                      << " was corrupted in the smoothed graph" << std::endl
                      << "original\t" << orig_seq << std::endl
                      << "smoothed\t" << smoothed_seq << std::endl;
            exit(1);
        }
        assert(orig_seq == smoothed_seq);
    });

    if (!consensus_mapping.empty()) {
        std::cerr << "[smoothxg::smooth_and_lace] sorting consensus"
                  << std::endl;
    }
    // consensus path and connections
    ips4o::parallel::sort(
        consensus_mapping.begin(), consensus_mapping.end(),
        [](const path_position_range_t &a, const path_position_range_t &b) {
            auto &a_id = as_integer(a.base_path);
            auto &b_id = as_integer(b.base_path);
            return (a_id < b_id || a_id == b_id && a.start_pos < b.start_pos);
        });

    // by definition, the consensus paths are embedded in our blocks, which
    // simplifies things we'll still need to add a new path for each consensus
    // path
    if (!consensus_mapping.empty()) {
        std::cerr << "[smoothxg::smooth_and_lace] embedding consensus"
                  << std::endl;
    }
    for (auto &pos_range : consensus_mapping) {
        auto &block = block_graphs[pos_range.target_graph_id];
        path_handle_t smoothed_path = smoothed.create_path_handle(
            block.get_path_name(pos_range.target_path));
        auto &id_trans = id_mapping[pos_range.target_graph_id];
        block.for_each_step_in_path(
            pos_range.target_path, [&](const step_handle_t &step) {
                handle_t h = block.get_handle_of_step(step);
                handle_t t = smoothed.get_handle(block.get_id(h) + id_trans,
                                                 block.get_is_reverse(h));
                smoothed.append_step(smoothed_path, t);
                // nb: by definition of our construction of smoothed
                // the consensus paths should have all their edges embedded
            });
    }
    // todo: validate the consensus paths as well
    return smoothed;
}

void write_gfa(std::unique_ptr<spoa::Graph> &graph, std::ostream &out,
               const std::vector<std::string> &sequence_names,
               bool include_consensus) {

    auto &nodes = graph->nodes();
    std::vector<std::int32_t> in_consensus(nodes.size(), -1);
    std::int32_t rank = 0;
    auto consensus = graph->consensus();
    for (const auto &id : consensus) {
        in_consensus[id] = rank++;
    }

    out << "H"
        << "\t"
        << "VN:Z:1.0" << std::endl;

    for (std::uint32_t i = 0; i < nodes.size(); ++i) {
        out << "S"
            << "\t" << i + 1 << "\t"
            << static_cast<char>(graph->decoder(nodes[i]->code()));
        if (in_consensus[i] != -1) {
            out << "\t"
                << "ic:Z:true";
        }
        out << std::endl;
        for (const auto &edge : nodes[i]->out_edges()) {
            out << "L"
                << "\t" << i + 1 << "\t"
                << "+"
                << "\t" << edge->end_node_id() + 1 << "\t"
                << "+"
                << "\t"
                << "0M"
                << "\t"
                << "ew:f:" << edge->total_weight();
            if (in_consensus[i] + 1 == in_consensus[edge->end_node_id()]) {
                out << "\t"
                    << "ic:Z:true";
            }
            out << std::endl;
        }
    }

    for (std::uint32_t i = 0; i < sequence_names.size(); ++i) {
        out << "P"
            << "\t" << sequence_names[i] << "\t";
        std::uint32_t node_id = graph->sequences_begin_nodes_ids()[i];
        while (true) {
            out << node_id + 1 << "+";
            if (!nodes[node_id]->successor(node_id, i)) {
                break;
            } else {
                out << ",";
            }
        }
        out << "\t"
            << "*" << std::endl;
    }

    if (include_consensus) {
        out << "P"
            << "\t"
            << "Consensus"
            << "\t";
        for (auto &id : graph->consensus()) {
            out << id + 1 << "+";
        }
        out << "\t"
            << "*" << std::endl;
    }
}

void build_odgi_abPOA(abpoa_t *ab, abpoa_para_t *abpt, odgi::graph_t &output,
                      const std::vector<std::string> &sequence_names,
                      const std::vector<bool> &aln_is_reverse,
                      const std::string &consensus_name,
                      bool include_consensus) {
    // std::cerr << "ENTERED build_odgi_abPOA" << std::endl;
    abpoa_graph_t *abg = ab->abg;
    if (abg->node_n <= 2) // how would this happen, and can we manage the error externally?
        return;
    int seq_n = sequence_names.size();

    // traverse graph
    int *in_degree = (int *)_err_malloc(abg->node_n * sizeof(int));
    int **read_paths = (int **)_err_malloc(seq_n * sizeof(int *)),
        *read_path_i = (int *)_err_calloc(seq_n, sizeof(int));
    int i, j, cur_id, pre_id, out_id, *id;
    for (i = 0; i < abg->node_n; ++i)
        in_degree[i] = abg->node[i].in_edge_n;
    for (i = 0; i < seq_n; ++i)
        read_paths[i] = (int *)_err_malloc(abg->node_n * sizeof(int));
    kdq_int_t *q = kdq_init_int();

    // Breadth-First-Search
    kdq_push_int(q, ABPOA_SRC_NODE_ID);
    while ((id = kdq_shift_int(q)) != 0) {
        cur_id = *id;
        if (cur_id == ABPOA_SINK_NODE_ID) {
            kdq_destroy_int(q);
            break;
        } else {
            if (cur_id != ABPOA_SRC_NODE_ID) {
                // output node
                // fprintf(stdout, "S\t%d\t%c\n", cur_id - 1,
                // "ACGTN"[abg->node[cur_id].base]); add node to output graph
                std::string seq = std::string(
                    1, static_cast<char>("ACGTN"[abg->node[cur_id].base]));
                // std::cerr << "seq: " << seq << std::endl;
                output.create_handle(seq, cur_id - 1);
                // std::cerr << "cur_id: " << (cur_id - 1) << std::endl;
                // output all links based pre_ids
                for (i = 0; i < abg->node[cur_id].in_edge_n; ++i) {
                    pre_id = abg->node[cur_id].in_id[i];
                    if (pre_id != ABPOA_SRC_NODE_ID) {
                        // output edge
                        // fprintf(stdout, "L\t%d\t+\t%d\t+\t0M\n", pre_id - 1,
                        // cur_id - 1); std::cerr << "cur_id edge: " << (cur_id
                        // - 1) << std::endl; std::cerr << "pre_id edge: " <<
                        // (pre_id - 1) << std::endl;
                        output.create_edge(output.get_handle(pre_id - 1),
                                           output.get_handle(cur_id - 1));
                    }
                }
                // add node id to read path
                int b, read_id;
                uint64_t num, tmp;
                b = 0;
                for (i = 0; i < abg->node[cur_id].read_ids_n; ++i) {
                    num = abg->node[cur_id].read_ids[i];
                    while (num) {
                        tmp = num & -num;
                        read_id = ilog2_64(abpt, tmp);
                        read_paths[b+read_id][read_path_i[b+read_id]++] = cur_id - 1;
                        num ^= tmp;
                    }
                    b += 64;
                }
            }
            for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
                out_id = abg->node[cur_id].out_id[i];
                if (--in_degree[out_id] == 0) {
                    kdq_push_int(q, out_id);
                }
            }
        }
    }

    // output read paths

    for (i = 0; i < seq_n; ++i) {
        path_handle_t p;
        std::vector<handle_t> steps;
        std::uint32_t node_id;
        if (!sequence_names.empty()) {
            // fprintf(stdout, "P\t%s\t", sequence_names[i]);
            // std::cerr << "P\t" << sequence_names[i] << "\t";
            p = output.create_path_handle(sequence_names[i]);
        } else {
            // fprintf(stdout, "P\t%d\t", i+1);
            p = output.create_path_handle(std::to_string(i + 1));
        }
        if (aln_is_reverse[i]) {
            for (j = read_path_i[i] - 1; j >= 0; --j) {
                // fprintf(stdout, "%d-", read_paths[i][j]);
                node_id = read_paths[i][j];
                output.append_step(p, output.flip(output.get_handle(node_id)));
                /*
                if (j != 0) {
                    fprintf(stdout, ",");
                } else {
                    fprintf(stdout, "\t*\n");
                }
                 */
            }
        } else {
            for (j = 0; j < read_path_i[i]; ++j) {
                // fprintf(stdout, "%d+", read_paths[i][j]);
                node_id = read_paths[i][j];
                output.append_step(p, output.get_handle(node_id));
                /*
                if (j != read_path_i[i]-1) {
                    fprintf(stdout, ",");
                }
                else {
                    fprintf(stdout, "\t*\n");
                }
                */
            }
        }
    }
    if (include_consensus) {
        // we already did that!
        // abpoa_generate_consensus(ab, abpt, seq_n, NULL, NULL, NULL, NULL,
        // NULL);
        int id = abg->node[ABPOA_SRC_NODE_ID].max_out_id;
        // fprintf(stdout, "P\tConsensus_sequence\t");
        path_handle_t p = output.create_path_handle(consensus_name);
        while (true) {
            // fprintf(stdout, "%d+", id-1);
            output.append_step(p, output.get_handle(id - 1));
            id = abg->node[id].max_out_id;
            if (id == ABPOA_SINK_NODE_ID) {
                break;
            }
            /*
            if (id != ABPOA_SINK_NODE_ID) {
                fprintf(stdout, ",");
            }
            else {
                // fprintf(stdout, "\t*\n");
                break;
            }
             */
        }
    }

    free(in_degree);
    for (i = 0; i < seq_n; ++i)
        free(read_paths[i]);
    free(read_paths);
    free(read_path_i);
    // std::cerr << "LEFT build_odgi_abPOA" << std::endl;
}

void build_odgi(std::unique_ptr<spoa::Graph> &graph, odgi::graph_t &output,
                const std::vector<std::string> &sequence_names,
                const std::vector<bool> &aln_is_reverse,
                const std::string &consensus_name, bool include_consensus) {

    auto &nodes = graph->nodes();
    std::vector<std::int32_t> in_consensus(nodes.size(), -1);
    std::int32_t rank = 0;
    auto consensus = graph->consensus();
    for (const auto &id : consensus) {
        in_consensus[id] = rank++;
    }

    for (std::uint32_t i = 0; i < nodes.size(); ++i) {
        std::string seq =
            std::string(1, static_cast<char>(graph->decoder(nodes[i]->code())));
        output.create_handle(seq, i + 1);
    }

    for (std::uint32_t i = 0; i < nodes.size(); ++i) {
        for (const auto &edge : nodes[i]->out_edges()) {
            output.create_edge(output.get_handle(i + 1),
                               output.get_handle(edge->end_node_id() + 1));
        }
    }

    for (std::uint32_t i = 0; i < sequence_names.size(); ++i) {
        path_handle_t p = output.create_path_handle(sequence_names[i]);
        std::uint32_t node_id = graph->sequences_begin_nodes_ids()[i];
        std::vector<handle_t> steps;
        while (true) {
            steps.push_back(output.get_handle(node_id + 1));
            if (!nodes[node_id]->successor(node_id, i)) {
                break;
            }
        }
        if (aln_is_reverse[i]) {
            for (auto handle_itr = steps.rbegin(); handle_itr != steps.rend();
                 ++handle_itr) {
                output.append_step(p, output.flip(*handle_itr));
            }
        } else {
            for (auto &handle : steps) {
                output.append_step(p, handle);
            }
        }
    }

    if (include_consensus) {
        path_handle_t p = output.create_path_handle(consensus_name);
        for (auto &id : graph->consensus()) {
            output.append_step(p, output.get_handle(id + 1));
        }
    }
}

} // namespace smoothxg
