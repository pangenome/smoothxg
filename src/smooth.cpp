#include "smooth.hpp"
#include "deps/abPOA/src/abpoa_align.c"
#include "deps/abPOA/src/abpoa_output.c"
#include "rkmh.hpp"

namespace smoothxg {

/*
void _clear_maf_block(ska::flat_hash_map<std::string, std::vector<maf_partial_row_t>> &maf){
    for (const auto &path_to_maf_rows : maf){
        for (auto &maf_row : path_to_maf_rows.second) {
            clear_string(maf_row.aligned_seq);
        }
    }

    maf.clear();
    ska::flat_hash_map<std::string, std::vector<maf_partial_row_t>>().swap(maf);
}

// if we want a vectorized layout representation of the block
void write_tsv_for_block(const xg::XG &graph,
                         const block_t &block,
                         const uint64_t &block_id,
                         const std::vector<std::string>& seqs,
                         const std::vector<std::string>& names) {
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
        ska::flat_hash_map<uint64_t, uint64_t> visits;
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
*/

void write_fasta_for_block(const xg::XG &graph,
                           const block_t &block,
                           const uint64_t &block_id,
                           const std::vector<std::string>& seqs,
                           const std::vector<std::string>& names,
                           const std::string& prefix,
                           const std::string& suffix) {
    std::string s = prefix + std::to_string(block_id) + suffix + ".fa";
    std::ofstream fasta(s.c_str());
    for (uint64_t i = 0; i < seqs.size(); ++i) {
        fasta << ">" << names[i] << " " << seqs[i].size() << std::endl
              << seqs[i] << std::endl;
    }
    fasta.close();
}

void append_to_sequence(const xg::XG &graph,
                        const path_handle_t &path_handle, const step_handle_t& starting_step,
                        std::basic_string<char> &seq, uint64_t &fwd_bp, uint64_t &rev_bp,
                        int poa_padding, bool on_the_left) {

    step_handle_t step = starting_step;
    const step_handle_t final_step = on_the_left ? graph.path_begin(path_handle) : graph.path_end(path_handle);

    uint64_t characters_to_add = poa_padding;
    std::string tmp;
    //ToDo: check if the condition is right
    while (step != final_step && characters_to_add > 0){
        const auto h = graph.get_handle_of_step(step);
        const auto l = graph.get_length(h);

        uint64_t characters_added;
        if (l <= characters_to_add) {
            // Take the full seq
            tmp.append(graph.get_sequence(h));
            characters_added = l;
        }else {
            // Take only the characters needed
            graph.get_sequence(h).substr(graph.get_sequence(h).size() - characters_to_add);
            tmp.append(
                    graph.get_sequence(h).substr(graph.get_sequence(h).size() - characters_to_add)
            );
            characters_added = characters_to_add;
        }
        if (graph.get_is_reverse(h)) {
            rev_bp += characters_added;
        } else {
            fwd_bp += characters_added;
        }
        characters_to_add -= characters_added;

        step = on_the_left ? graph.get_previous_step(step) : graph.get_next_step(step);
    }

    if (on_the_left){
        //std::cerr << graph.get_path_name(path_handle) << " - HEAD: ";
        while (characters_to_add > 0) {
            seq.append("N");
            //std::cerr << "N";
            --characters_to_add;
        }
        //std::cerr << tmp << std::endl;
        seq.append(tmp);
    } else {
        //std::cerr << graph.get_path_name(path_handle) << " - TAIL: " << tmp;
        seq.append(tmp);

        while (characters_to_add > 0) {
            seq.append("N");
            //std::cerr << "N";
            --characters_to_add;
        }
        //std::cerr << "\n";
    }
}

odgi::graph_t* smooth_abpoa(const xg::XG &graph, const block_t &block, const uint64_t block_id,
                            int poa_m, int poa_n, int poa_g,
                            int poa_e, int poa_q, int poa_c,
                            int poa_padding,
                            bool local_alignment,
                            std::unique_ptr<ska::flat_hash_map<std::string, std::vector<maf_partial_row_t>>>& maf, bool keep_sequence,
                            bool banded_alignment,
                            const std::string &consensus_name,
                            bool save_block_fastas) {
    // collect sequences
    std::vector<std::string> seqs;
    std::vector<std::string> names;
    std::vector<bool> is_rev;
    std::size_t max_sequence_size = 0;

    for (auto &path_range : block.path_ranges) {
        seqs.emplace_back();
        auto &seq = seqs.back();
        uint64_t fwd_bp = 0;
        uint64_t rev_bp = 0;
        const path_handle_t path_handle = graph.get_path_handle_of_step(path_range.begin);

        append_to_sequence(graph,
                           path_handle, path_range.begin,
                           seq, fwd_bp, rev_bp,
                           poa_padding, true);

        for (step_handle_t step = path_range.begin; step != path_range.end;
             step = graph.get_next_step(step)) {
            const auto h = graph.get_handle_of_step(step);
            const auto l = graph.get_length(h);
            seq.append(graph.get_sequence(h));
            if (graph.get_is_reverse(h)) {
                rev_bp += l;
            } else {
                fwd_bp += l;
            }
        }

        append_to_sequence(graph,
                           path_handle, path_range.end,
                           seq, fwd_bp, rev_bp,
                           poa_padding, false);

        is_rev.push_back(rev_bp > fwd_bp);
        if (is_rev.back()) {
            odgi::reverse_complement_in_place(seqs.back());
        }
        std::stringstream namess;
        namess << graph.get_path_name(path_handle)
               << "_" << graph.get_position_of_step(path_range.begin);
        names.push_back(namess.str());

        max_sequence_size = std::max(max_sequence_size, seq.size());
    }

    auto start_time = std::chrono::steady_clock::now();

    if (save_block_fastas) {
        write_fasta_for_block(graph, block, block_id, seqs, names, "smoothxg_into_abpoa_pad" + std::to_string(poa_padding) + "_");
    }

    auto* output_graph = new odgi::graph_t();

    // if the graph would be empty, bail out
    if (max_sequence_size == 0) {
        return output_graph;
    }

    const bool add_consensus = !consensus_name.empty();

    // initialize abPOA parameters
    abpoa_para_t *abpt = abpoa_init_para();

    // if we want to do local alignments
    if (local_alignment) {
        abpt->align_mode = ABPOA_LOCAL_MODE;
    } else {
        abpt->align_mode = ABPOA_GLOBAL_MODE;
    }
    // abpt->zdrop = -1;     // disable zdrop
    // abpt->end_bonus = -1; // disable end bouns
    if (!banded_alignment) {
        abpt->wb = -1;
    } else {
        abpt->wb = 311;
    }
    abpt->wf = 0.03; // hmm
    abpt->amb_strand = 0; //we align based on orientation relative to the graph
    // abpt->ret_cigar = 1;  // return cigar
    abpt->rev_cigar = 0;  // reverse cigar
    abpt->out_cons = add_consensus;
    // abpt->out_fq = 0;     // output consensus sequence in fastq
    abpt->out_gfa = 1; // must be set to get the graph
    abpt->out_msa = maf != nullptr ? 1 : 0; // must be set when we extract the MSA
    // abpt->max_n_cons = 1; // number of max. generated consensus sequence

    // score matrix
    abpt->match = poa_m;
    abpt->mismatch = poa_n;
    abpt->gap_open1 = poa_g;
    abpt->gap_open2 = poa_q;
    abpt->gap_ext1 = poa_e;
    abpt->gap_ext2 = poa_c;

    abpt->disable_seeding = 1; // disable seeding (allow seeding greatly reduces runtime and memory, but reducing the accuracy)

    abpt->k = 19;
    abpt->w = 10;
    abpt->min_w = 3313;
    // abpt->progressive_poa = 0; // progressive partial order alignment

    // finalize parameters
    abpoa_post_set_para(abpt);

    // initialize abPOA
    abpoa_t *ab = abpoa_init();

    // alloc our seq count
     // collect sequence length, transform ACGT to 0123
    const int n_seq = seqs.size();
    int *seq_lens = (int *)_err_malloc(sizeof(int) * n_seq);
    auto **bseqs = (uint8_t **)_err_malloc(sizeof(uint8_t *) * n_seq);
    for (int i = 0; i < n_seq; ++i) {
        seq_lens[i] = seqs[i].size();
        bseqs[i] = (uint8_t *)_err_malloc(sizeof(uint8_t) * seq_lens[i]);
        for (int j = 0; j < seq_lens[i]; ++j) {
            bseqs[i][j] = ab_char26_table[(int)seqs[i][j]];
        }
    }

    abpoa_reset(ab, abpt, 1024);
    abpoa_seq_t *abs = ab->abs; int i, exist_n_seq = ab->abs->n_seq;

    // alloc our seq count
    abs->n_seq += n_seq; abpoa_realloc_seq(abs);

    // no seq names
    for (i = 0; i < n_seq; ++i) {
        abs->name[exist_n_seq+i].l = 0; abs->name[exist_n_seq+i].m = 0;
    }

    // determine the max sequence length
    // int max_len = 0;
    // for (i = 0; i < n_seq; ++i) {
    //     if (seq_lens[i] > max_len) max_len = seq_lens[i];
    // }

    abpoa_poa(ab, abpt, bseqs, seq_lens, exist_n_seq, n_seq);

    // It seems not necessary, as the Breadth-First-Search (in build_odgi_abPOA) follows the graph topology
    //abpoa_topological_sort(ab->abg, abpt);

    if (maf != nullptr){
        abpoa_generate_rc_msa(ab, abpt);
    }

    if (add_consensus) {
        abpoa_generate_consensus(ab, abpt);
        if (ab->abg->is_called_cons == 0) {
            err_printf("ERROR: no consensus sequence generated.\n");
            exit(1);
        }
        is_rev.push_back(false);  // the consensus is considered in forward

        //abpoa_output_fx_consensus(ab, abpt, stdout);
    }

    if (maf != nullptr) {
        //abpoa_output_rc_msa(ab, abpt, stdout);

        for (int i = 0; i < n_seq; ++i) {
            // Remove padded characters
            int j = 0;
            uint64_t characters_to_remove = poa_padding;
            while (characters_to_remove > 0) {
                if (ab->abc->msa_base[i][j] != 5){
                    ab->abc->msa_base[i][j] = 5;
                    --characters_to_remove;
                }

                ++j;
            }

            j = ab->abc->msa_len;
            characters_to_remove = poa_padding;
            while (characters_to_remove > 0) {
                --j;
                if (ab->abc->msa_base[i][j] != 5){
                    ab->abc->msa_base[i][j] = 5;
                    --characters_to_remove;
                }
            }
        }

        // Find the starting position where to trim the MSA
        uint64_t start_pos_to_trim = 0;
        for (; start_pos_to_trim < ab->abc->msa_len; ++start_pos_to_trim) {
            int i = 0;
            // Find a non-gap character in the current MSA-column
            for (; i < n_seq; ++i) {
                if (ab->abc->msa_base[i][start_pos_to_trim] != 5){
                    break;
                }
            }

            if (i < n_seq) {
                // A non-gap character was found
                break;
            }
        }

        // Find the ending position where to trim the MSA
        int64_t end_pos_to_trim = ab->abc->msa_len - 1;
        for (; end_pos_to_trim >= 0; --end_pos_to_trim) {
            int i = 0;
            // Find a non-gap character in the current MSA-column
            for (; i < n_seq; ++i) {
                if (ab->abc->msa_base[i][end_pos_to_trim] != 5){
                    break;
                }
            }

            if (i < n_seq) {
                // A non-gap character was found
                break;
            }
        }
        end_pos_to_trim += 1;


        uint64_t num_seqs = n_seq + (add_consensus ? 1 : 0);

        for(uint64_t seq_rank = 0; seq_rank < num_seqs; seq_rank++) {
            std::basic_string<char> aligned_seq;

            std::string path_name;
            uint64_t seq_size;
            uint64_t path_length;
            uint64_t record_start;
            if (!add_consensus || seq_rank < num_seqs - 1) {
                if (keep_sequence){
                    for (int j = start_pos_to_trim; j < end_pos_to_trim; ++j) {
                        aligned_seq += ab_char256_table[ab->abc->msa_base[seq_rank][j]];
                    }
                }

                auto path_handle = graph.get_path_handle_of_step(block.path_ranges[seq_rank].begin);

                path_name = graph.get_path_name(path_handle);
                path_length = graph.get_path_length(path_handle);

                // If the strand field is "-" then this is the start relative to the reverse-complemented source sequence
                uint64_t path_range_begin = graph.get_position_of_step(block.path_ranges[seq_rank].begin);
                auto last_step = graph.get_previous_step(block.path_ranges[seq_rank].end);
                record_start = is_rev[seq_rank] ?
                               (path_length - graph.get_position_of_step(last_step) -
                                graph.get_length(graph.get_handle_of_step(last_step))) :
                               path_range_begin;

                seq_size = seqs[seq_rank].size() - 2 * poa_padding; // <==> block.path_ranges[seq_rank].length
            } else {
                // The last sequence is the gapped consensus
                if (keep_sequence){
                    // TODO: is it really enough to replace the old code?
                    for (i = 0; i < ab->abc->msa_len; ++i){
                        aligned_seq += ab_char256_table[ab->abc->msa_base[ab->abc->n_seq][i]];
                    }

                    // int j, k, aligned_id, rank;
                    // i = ab->abg->node[ABPOA_SRC_NODE_ID].max_out_id;
                    // int last_rank = 1;
                    // while (i != ABPOA_SINK_NODE_ID) {
                    //     rank = abpoa_graph_node_id_to_msa_rank(ab->abg, i);
                    //     for (k = 0; k < ab->abg->node[i].aligned_node_n; ++k) {
                    //         aligned_id = ab->abg->node[i].aligned_node_id[k];
                    //         rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(ab->abg, aligned_id));
                    //     }
                    //     // last_rank -> rank : -
                    //     for (k = last_rank; k < rank; ++k) aligned_seq += '-';
                    //     // rank : base
                    //     aligned_seq += "ACGTN"[ab->abg->node[i].base];
                    //     last_rank = rank+1;
                    //     i = ab->abg->node[i].ma.max_out_id;
                    // }
                    // // last_rank -> msa_l:-
                    // for (k = last_rank; k <= ab->abc->msa_len; ++k) aligned_seq += '-';

                    aligned_seq = aligned_seq.substr(start_pos_to_trim, end_pos_to_trim - start_pos_to_trim);
                }

                path_name = consensus_name;
                path_length = ab->abc->cons_len[0] - 2 * poa_padding;
                record_start = 0;
                seq_size = path_length;
            }

            if (maf->count(path_name) == 0) {
                maf->insert(std::pair<std::string, std::vector<maf_partial_row_t>>(
                        path_name,
                        std::vector<maf_partial_row_t>()
                ));
            }

            (*maf)[path_name].push_back({
                record_start,
                seq_size,
                is_rev[seq_rank] == 1,
                path_length,
                aligned_seq
            });

            clear_string(path_name);
            clear_string(aligned_seq);
        }
    }

    // free memory
    for (i = 0; i < n_seq; ++i) {
        free(bseqs[i]);
    }
    free(bseqs);
    free(seq_lens);

    if (save_block_fastas) {
        std::chrono::duration<double> elapsed_time = std::chrono::steady_clock::now() - start_time;
        write_fasta_for_block(graph, block, block_id, seqs, names, "smoothxg_into_abpoa_pad" + std::to_string(poa_padding) + "_", "_in_" +  std::to_string(elapsed_time.count()) + "s");
    }

    odgi::graph_t block_graph;
    build_odgi_abPOA(ab, abpt, &block_graph, names, is_rev, consensus_name, poa_padding, add_consensus);
    abpoa_free(ab);
    abpoa_free_para(abpt);

    // normalize the representation, allowing for nodes > 1bp
    //auto graph_copy = output_graph;
    if (!odgi::algorithms::unchop(block_graph)) {
        std::cerr << "[smoothxg::smooth_abpoa] error: unchop failure, saving before/after graph to disk" << std::endl;
        //std::ofstream a("smoothxg_unchop_failure_before.gfa");
        //graph_copy.to_gfa(a);
        //a.close();
        std::ofstream b("smoothxg_unchop_failure_after.gfa");
        block_graph.to_gfa(b);
        b.close();
        exit(1);
    }

    // order the graph
    block_graph.apply_ordering(odgi::algorithms::topological_order(&block_graph), true);

    // copy the graph to avoid memory fragmentation issues
    block_graph.for_each_handle(
            [&](const handle_t& old_handle) {
                output_graph->create_handle(
                        block_graph.get_sequence(old_handle),
                        block_graph.get_id(old_handle));
            });

    block_graph.for_each_handle(
            [&](const handle_t& curr) {
                block_graph.follow_edges(
                        curr, false,
                        [&](const handle_t& next) {
                            output_graph->create_edge(
                                    output_graph->get_handle(block_graph.get_id(curr),
                                                             block_graph.get_is_reverse(curr)),
                                    output_graph->get_handle(block_graph.get_id(next),
                                                             block_graph.get_is_reverse(next)));
                        });
                block_graph.follow_edges(
                        curr, true,
                        [&](const handle_t& prev) {
                            output_graph->create_edge(
                                    output_graph->get_handle(block_graph.get_id(prev),
                                                             block_graph.get_is_reverse(prev)),
                                    output_graph->get_handle(block_graph.get_id(curr),
                                                             block_graph.get_is_reverse(curr)));
                        });
            });

    block_graph.for_each_path_handle(
            [&](const path_handle_t& old_path) {
                path_handle_t new_path = output_graph->create_path_handle(block_graph.get_path_name(old_path));
                block_graph.for_each_step_in_path(old_path, [&](const step_handle_t& step) {
                    handle_t old_handle = block_graph.get_handle_of_step(step);
                    handle_t new_handle = output_graph->get_handle(
                            block_graph.get_id(old_handle),
                            block_graph.get_is_reverse(old_handle));
                    output_graph->append_step(new_path, new_handle);
                });
            });

    return output_graph;
}

odgi::graph_t* smooth_spoa(const xg::XG &graph, const block_t &block,
                           const uint64_t block_id,
                           std::int8_t poa_m, std::int8_t poa_n, std::int8_t poa_g,
                           std::int8_t poa_e, std::int8_t poa_q, std::int8_t poa_c,
                           int poa_padding,
                           bool local_alignment,
                           std::unique_ptr<ska::flat_hash_map<std::string, std::vector<maf_partial_row_t>>>& maf, bool keep_sequence,
                           const std::string &consensus_name,
                           bool save_block_fastas) {
    // collect sequences
    std::vector<std::string> seqs;
    std::vector<std::string> names;
    std::vector<bool> is_rev;
    std::size_t max_sequence_size = 0;

    for (auto &path_range : block.path_ranges) {
        seqs.emplace_back();
        auto &seq = seqs.back();
        uint64_t fwd_bp = 0;
        uint64_t rev_bp = 0;
        const path_handle_t path_handle = graph.get_path_handle_of_step(path_range.begin);

        append_to_sequence(graph,
                           path_handle, path_range.begin,
                           seq, fwd_bp, rev_bp,
                           poa_padding, true);

        for (step_handle_t step = path_range.begin; step != path_range.end;
             step = graph.get_next_step(step)) {
            const auto h = graph.get_handle_of_step(step);
            const auto l = graph.get_length(h);
            seq.append(graph.get_sequence(h));
            if (graph.get_is_reverse(h)) {
                rev_bp += l;
            } else {
                fwd_bp += l;
            }
        }

        append_to_sequence(graph,
                           path_handle, path_range.end,
                           seq, fwd_bp, rev_bp,
                           poa_padding, false);

        is_rev.push_back(rev_bp > fwd_bp);
        if (is_rev.back()) {
            odgi::reverse_complement_in_place(seqs.back());
        }
        std::stringstream namess;
        namess << graph.get_path_name(path_handle)
               << "_" << graph.get_position_of_step(path_range.begin);
        names.push_back(namess.str());

        max_sequence_size = std::max(max_sequence_size, seq.size());
    }

    if (save_block_fastas) {
        write_fasta_for_block(graph, block, block_id, seqs, names, "smoothxg_into_spoa");
    }

    auto* output_graph = new odgi::graph_t();

    // if the graph would be empty, bail out
    if (max_sequence_size == 0) {
        return output_graph;
    }

    auto spoa_algorithm = local_alignment ? spoa::AlignmentType::kSW : spoa::AlignmentType::kNW;
    auto alignment_engine = spoa::AlignmentEngine::Create(
        spoa_algorithm,
        poa_m, poa_n, poa_g, poa_e, poa_q, poa_c);

    spoa::Graph poa_graph{};

    int i = 0;
    for (auto &seq : seqs) {
        auto alignment = alignment_engine->Align(seq, poa_graph);
        try {
            // could give weight here to influence consensus
            poa_graph.AddAlignment(alignment, seq);
        } catch (std::invalid_argument &exception) {
            std::cerr << exception.what() << std::endl;
            assert(false);
        }
    }

    std::string consensus;
    if (!consensus_name.empty()){
        consensus = poa_graph.GenerateConsensus();
        is_rev.push_back(false);  // the consensus is considered in forward
    }

    if (maf != nullptr) {
        bool add_consensus = !consensus_name.empty();

        std::vector<std::string> msa =
            poa_graph.GenerateMultipleSequenceAlignment(add_consensus);

        const int msa_l = msa[0].size();

        for (int i = 0; i < msa.size(); ++i) {
            // Remove padded characters
            int j = 0;
            uint64_t characters_to_remove = poa_padding;
            while (characters_to_remove > 0) {
                if (msa[i][j] != 5){
                    msa[i][j] = 5;
                    --characters_to_remove;
                }

                ++j;
            }

            j = msa_l;
            characters_to_remove = poa_padding;
            while (characters_to_remove > 0) {
                --j;
                if (msa[i][j] != 5){
                    msa[i][j] = 5;
                    --characters_to_remove;
                }
            }
        }

        // Find the starting position where to trim the MSA
        uint64_t start_pos_to_trim = 0;
        for (; start_pos_to_trim < msa_l; ++start_pos_to_trim) {
            int i = 0;
            // Find a non-gap character in the current MSA-column
            for (; i < msa.size(); ++i) {
                if (msa[i][start_pos_to_trim] != 5){
                    break;
                }
            }

            if (i < msa.size()) {
                // A non-gap character was found
                break;
            }
        }

        // Find the ending position where to trim the MSA
        int64_t end_pos_to_trim = msa_l - 1;
        for (; end_pos_to_trim >= 0; --end_pos_to_trim) {
            int i = 0;
            // Find a non-gap character in the current MSA-column
            for (; i < msa.size(); ++i) {
                if (msa[i][end_pos_to_trim] != 5){
                    break;
                }
            }

            if (i < msa.size()) {
                // A non-gap character was found
                break;
            }
        }
        end_pos_to_trim += 1;


        uint64_t num_seqs = msa.size();
        for(uint64_t seq_rank = 0; seq_rank < num_seqs; seq_rank++){
            std::string path_name;
            uint64_t seq_size;
            uint64_t path_length;

            uint64_t record_start;
            if (!add_consensus || seq_rank < num_seqs - 1){
                auto path_handle = graph.get_path_handle_of_step(block.path_ranges[seq_rank].begin);

                path_name = graph.get_path_name(path_handle);
                path_length = graph.get_path_length(path_handle);

                // If the strand field is "-" then this is the start relative to the reverse-complemented source sequence
                uint64_t path_range_begin = graph.get_position_of_step(block.path_ranges[seq_rank].begin);
                auto last_step = graph.get_previous_step(block.path_ranges[seq_rank].end);
                record_start = is_rev[seq_rank] ?
                               (path_length - graph.get_position_of_step(last_step) - graph.get_length(graph.get_handle_of_step(last_step))):
                               path_range_begin;

                seq_size = seqs[seq_rank].size() - 2 * poa_padding;// <==> block.path_ranges[seq_rank].length
            }else{
                // The last sequence is the gapped consensus

                path_name = consensus_name;
                path_length = consensus.size() - 2 * poa_padding;
                record_start = 0;
                seq_size = path_length;
            }

            if (maf->count(path_name) == 0) {
                maf->insert(std::pair<std::string, std::vector<maf_partial_row_t>>(
                        path_name,
                        std::vector<maf_partial_row_t>()
                ));
            }

            (*maf)[path_name].push_back({
                    record_start,
                    seq_size,
                    is_rev[seq_rank],
                    path_length,
                    keep_sequence ? msa[seq_rank].substr(start_pos_to_trim, end_pos_to_trim - start_pos_to_trim) : ""
            });

            clear_string(path_name);
            clear_string(msa[seq_rank]);
        }

        clear_vector(msa);
    }

    // write the graph, with consensus as a path
    // odgi::graph_t output_graph;
    // convert the poa graph into our output format
    // poa_graph->print_gfa(std::cout, names, true);
    odgi::graph_t block_graph;
    build_odgi_SPOA(poa_graph, &block_graph, names, poa_padding, is_rev, consensus_name, !consensus_name.empty());

        // normalize the representation, allowing for nodes > 1bp
    //auto graph_copy = output_graph;
    if (!odgi::algorithms::unchop(block_graph)) {
        std::cerr << "[smoothxg::smooth_abpoa] error: unchop failure, saving before/after graph to disk" << std::endl;
        //std::ofstream a("smoothxg_unchop_failure_before.gfa");
        //graph_copy.to_gfa(a);
        //a.close();
        std::ofstream b("smoothxg_unchop_failure_after.gfa");
        block_graph.to_gfa(b);
        b.close();
        exit(1);
    }

    // order the graph
    block_graph.apply_ordering(odgi::algorithms::topological_order(&block_graph), true);

    // copy the graph to avoid memory fragmentation issues
    block_graph.for_each_handle(
            [&](const handle_t& old_handle) {
                output_graph->create_handle(
                        block_graph.get_sequence(old_handle),
                        block_graph.get_id(old_handle));
            });

    block_graph.for_each_handle(
            [&](const handle_t& curr) {
                block_graph.follow_edges(
                        curr, false,
                        [&](const handle_t& next) {
                            output_graph->create_edge(
                                    output_graph->get_handle(block_graph.get_id(curr),
                                                             block_graph.get_is_reverse(curr)),
                                    output_graph->get_handle(block_graph.get_id(next),
                                                             block_graph.get_is_reverse(next)));
                        });
                block_graph.follow_edges(
                        curr, true,
                        [&](const handle_t& prev) {
                            output_graph->create_edge(
                                    output_graph->get_handle(block_graph.get_id(prev),
                                                             block_graph.get_is_reverse(prev)),
                                    output_graph->get_handle(block_graph.get_id(curr),
                                                             block_graph.get_is_reverse(curr)));
                        });
            });

    block_graph.for_each_path_handle(
            [&](const path_handle_t& old_path) {
                path_handle_t new_path = output_graph->create_path_handle(block_graph.get_path_name(old_path));
                block_graph.for_each_step_in_path(old_path, [&](const step_handle_t& step) {
                    handle_t old_handle = block_graph.get_handle_of_step(step);
                    handle_t new_handle = output_graph->get_handle(
                            block_graph.get_id(old_handle),
                            block_graph.get_is_reverse(old_handle));
                    output_graph->append_step(new_path, new_handle);
                });
            });

    // output_graph.to_gfa(out);
    return output_graph;
}

void _put_block_in_group(
        maf_t &merged_maf_blocks, uint64_t block_id, uint64_t num_seq_in_block,
        std::string consensus_name,
        std::vector<std::unique_ptr<ska::flat_hash_map<std::string, std::vector<maf_partial_row_t>>>> &mafs,
        bool new_block_on_the_left,
        bool flip_block_before_merging
){
    uint64_t alignment_size_merged_maf_blocks =
            merged_maf_blocks.block_ids.empty() ? 0 : merged_maf_blocks.rows.begin()->second.begin() ->aligned_seq.length();

    // If in the merged group there are sequences from previous blocks, put gaps for the new added paths from new blocks
    std::string gaps = "";
    for (uint64_t j = 0; j < alignment_size_merged_maf_blocks; j++){ gaps += "-"; }

    for (const auto& path_to_maf_rows : *mafs[block_id]) {
        // Do not check the consensus (always forward)
        if (path_to_maf_rows.first != consensus_name){
            if (merged_maf_blocks.rows.count(path_to_maf_rows.first) == 0) {
                merged_maf_blocks.rows.insert(std::pair<std::string, std::vector<maf_partial_row_t>>(
                        path_to_maf_rows.first,
                        std::vector<maf_partial_row_t>()
                ));

                for (auto& maf_row : path_to_maf_rows.second) {
                    uint64_t maf_row_record_start = flip_block_before_merging ? maf_row.path_length - (maf_row.record_start + maf_row.seq_size) : maf_row.record_start;

                    if (flip_block_before_merging) {
                        odgi::reverse_complement_in_place((std::string&)maf_row.aligned_seq);
                    }

                    merged_maf_blocks.rows[path_to_maf_rows.first].push_back({
                            maf_row_record_start,
                            maf_row.seq_size,
                            (bool)(flip_block_before_merging ^ maf_row.is_reversed),
                            maf_row.path_length,

                            new_block_on_the_left ? maf_row.aligned_seq + gaps : gaps + maf_row.aligned_seq
                    });
                }
            } else {
                // Try to merge all the mergeable intervals

                std::vector<uint64_t> unmerged_rows;

                for (uint64_t rank_row = 0; rank_row < path_to_maf_rows.second.size(); ++rank_row) {
                    //std::cerr << "rank_row " << rank_row << std::endl;

                    auto& maf_row = path_to_maf_rows.second[rank_row];
                    uint64_t maf_row_record_start = flip_block_before_merging ? maf_row.path_length - (maf_row.record_start + maf_row.seq_size) : maf_row.record_start;

                    bool merged = false;

                    for (auto &merged_maf_prow : merged_maf_blocks.rows[path_to_maf_rows.first]) {
                        // Check the length to avoid merging more rows from the same last block
                        if ((flip_block_before_merging ^ maf_row.is_reversed) == merged_maf_prow.is_reversed &&
                        merged_maf_prow.aligned_seq.length() == alignment_size_merged_maf_blocks){
                            if (merged_maf_prow.is_reversed){
                                if ((merged_maf_prow.path_length - merged_maf_prow.record_start) == (maf_row.path_length - (maf_row_record_start + maf_row.seq_size))) {
                                    // merged_maf_row_end == maf_row_begin, new row on the left

                                    merged_maf_prow.record_start -= maf_row.seq_size;

                                    if (flip_block_before_merging) {
                                        odgi::reverse_complement_in_place((std::string&)maf_row.aligned_seq);
                                    }

                                    merged_maf_prow.aligned_seq = maf_row.aligned_seq + merged_maf_prow.aligned_seq;
                                    merged_maf_prow.seq_size += maf_row.seq_size;

                                    merged = true;
                                    break;
                                } else if ((maf_row.path_length - maf_row_record_start) == (merged_maf_prow.path_length - (merged_maf_prow.record_start + merged_maf_prow.seq_size))) {
                                    // maf_row_end == merged_maf_row_begin, new row on the right

                                    if (flip_block_before_merging) {
                                        odgi::reverse_complement_in_place((std::string&)maf_row.aligned_seq);
                                    }

                                    merged_maf_prow.aligned_seq += maf_row.aligned_seq;
                                    merged_maf_prow.seq_size += maf_row.seq_size;

                                    merged = true;
                                    break;
                                }
                            } else {
                                if ((merged_maf_prow.record_start + merged_maf_prow.seq_size) == maf_row_record_start) {
                                    // merged_maf_row_end == maf_row_begin, new row on the right

                                    if (flip_block_before_merging) {
                                        odgi::reverse_complement_in_place((std::string&)maf_row.aligned_seq);
                                    }

                                    merged_maf_prow.aligned_seq += maf_row.aligned_seq;
                                    merged_maf_prow.seq_size += maf_row.seq_size;

                                    merged = true;
                                    break;
                                } else if ((maf_row_record_start + maf_row.seq_size) == merged_maf_prow.record_start) {
                                    // maf_row_end == merged_maf_row_begin, new row on the left

                                    merged_maf_prow.record_start -= maf_row.seq_size;

                                    if (flip_block_before_merging) {
                                        odgi::reverse_complement_in_place((std::string&)maf_row.aligned_seq);
                                    }

                                    merged_maf_prow.aligned_seq = maf_row.aligned_seq + merged_maf_prow.aligned_seq ;
                                    merged_maf_prow.seq_size += maf_row.seq_size;

                                    merged = true;
                                    break;
                                }
                            }
                        }
                    }

                    if (!merged) {
                        unmerged_rows.push_back(rank_row);
                    }
                }

                for (auto& rank_row : unmerged_rows) {
                    auto& maf_row = path_to_maf_rows.second[rank_row];
                    uint64_t maf_row_record_start = flip_block_before_merging ? maf_row.path_length - (maf_row.record_start + maf_row.seq_size) : maf_row.record_start;

                    if (flip_block_before_merging) {
                        odgi::reverse_complement_in_place((std::string&)maf_row.aligned_seq);
                    }

                    merged_maf_blocks.rows[path_to_maf_rows.first].push_back({
                        maf_row_record_start,
                        maf_row.seq_size,
                        (bool)(flip_block_before_merging ^ maf_row.is_reversed),
                        maf_row.path_length,
                        new_block_on_the_left ? maf_row.aligned_seq + gaps : gaps + maf_row.aligned_seq
                    });
                }
            }
        }
    }

    clear_string(gaps);

    // The merged consensus is created when the merged block is written into a file
    if (!consensus_name.empty()){
        // IMPORTANT: it assumes a single consensus sequence!
        auto &maf_row = (*mafs[block_id])[consensus_name][0];

        //uint64_t maf_row_record_start = flip_block_before_merging ? maf_row.path_length - (maf_row.record_start + maf_row.seq_size) : maf_row.record_start;
        if (flip_block_before_merging) {
            odgi::reverse_complement_in_place(maf_row.aligned_seq);
        }

        if (new_block_on_the_left){
            merged_maf_blocks.consensus_rows.push_front(std::pair<std::string, maf_partial_row_t>(
                    consensus_name,
                    {
                            maf_row.record_start,//maf_row_record_start,
                            maf_row.seq_size,
                            maf_row.is_reversed,//(bool) (flip_block_before_merging ^ maf_row.is_reversed),
                            maf_row.path_length,
                            maf_row.aligned_seq
                    }
            ));
        }else {
            merged_maf_blocks.consensus_rows.push_back(std::pair<std::string, maf_partial_row_t>(
                    consensus_name,
                    {
                            maf_row.record_start,//maf_row_record_start,
                            maf_row.seq_size,
                            maf_row.is_reversed,//(bool) (flip_block_before_merging ^ maf_row.is_reversed),
                            maf_row.path_length,
                            maf_row.aligned_seq
                    }
            ));
        }
    }

    // Put gaps for paths not present in the last merged block (block_id) respect to the merged group

    // I take the length from a one of the path present in the last merged block
    uint64_t num_gaps_to_add = mafs[block_id]->begin()->second[0].aligned_seq.size();
    alignment_size_merged_maf_blocks += num_gaps_to_add;

    for (uint64_t  i = 0; i < num_gaps_to_add; i++){ gaps += "-"; }

    for (auto &path_to_maf_rows_m : merged_maf_blocks.rows){
        for (auto &merged_maf_prow : path_to_maf_rows_m.second){
            if (merged_maf_prow.aligned_seq.length() < alignment_size_merged_maf_blocks){
                if (new_block_on_the_left){
                    merged_maf_prow.aligned_seq = gaps + merged_maf_prow.aligned_seq;
                } else {
                    merged_maf_prow.aligned_seq += gaps;
                }

                /*if (merged_maf_prow.is_reversed){
                    if (new_block_on_the_left){
                        merged_maf_prow.aligned_seq += gaps;
                    } else {
                        merged_maf_prow.aligned_seq = gaps + merged_maf_prow.aligned_seq;
                    }
                }else{
                    if (new_block_on_the_left){
                        merged_maf_prow.aligned_seq = gaps + merged_maf_prow.aligned_seq;
                    } else {
                        merged_maf_prow.aligned_seq += gaps;
                    }
                }*/
            }
        }
    }

    clear_string(gaps);

    if (new_block_on_the_left){
        merged_maf_blocks.block_ids.insert(merged_maf_blocks.block_ids.begin(), block_id);
    } else {
        merged_maf_blocks.block_ids.push_back(block_id);
    }
}

void _write_merged_maf_blocks(
        maf_t& merged_maf_blocks,
        ska::flat_hash_set<uint64_t> &inverted_merged_block_id_intervals_ranks,
        std::vector<IITree<uint64_t, uint64_t>> &merged_block_id_intervals_tree_vector,
        std::vector<std::string> &block_id_ranges_vector,
        std::vector<bool> &is_block_in_a_merged_group,
        bool produce_maf, std::ofstream& out_maf,
        bool add_consensus, std::string consensus_base_name,
        bool fraction_below_threshold,
        bool preserve_unmerged_consensus
        ) {
    uint64_t merged_maf_blocks_size = merged_maf_blocks.block_ids.size();

    /*std::stringstream joined_block_ids;
    for (size_t i = 0; i < merged_maf_blocks_size; ++i) {
        if (i != 0) { joined_block_ids << ","; }
        joined_block_ids << merged_maf_blocks.block_ids[i];
    }*/

    // get the min/max
    uint64_t min_block_id = merged_maf_blocks.block_ids.front();
    uint64_t max_block_id = merged_maf_blocks.block_ids.back();
    if (min_block_id > max_block_id) {
        // It means that the blocks have been joined from the left
        min_block_id = max_block_id;
        max_block_id = merged_maf_blocks.block_ids.front();
    }

    std::string block_id_ranges = std::to_string(min_block_id);
    std::string full_block_id_ranges = std::to_string(merged_maf_blocks.block_ids.front());
    if (merged_maf_blocks_size > 1) {
        full_block_id_ranges = "";
        block_id_ranges += "-" + std::to_string(max_block_id);

        bool inverted_merged_block_ids = merged_maf_blocks.block_ids.front() > merged_maf_blocks.block_ids.back();

        if (inverted_merged_block_ids) {
            inverted_merged_block_id_intervals_ranks.insert(merged_block_id_intervals_tree_vector.size());
        }

        merged_block_id_intervals_tree_vector.emplace_back();
        auto& merged_block_id_intervals_tree = merged_block_id_intervals_tree_vector.back();

        if (add_consensus) {
            is_block_in_a_merged_group[merged_maf_blocks.block_ids[0]] = true;
        }

        uint64_t begin_pos = 0;
        for (uint64_t i = 1; i < merged_maf_blocks_size; ++i) {
            bool contiguos = inverted_merged_block_ids ?
                             merged_maf_blocks.block_ids[i - 1] - merged_maf_blocks.block_ids[i] == 1 :
                             merged_maf_blocks.block_ids[i] - merged_maf_blocks.block_ids[i - 1] == 1;

            if (!contiguos) {
                if (inverted_merged_block_ids){
                    merged_block_id_intervals_tree.add(
                            merged_maf_blocks.block_ids[i - 1],
                            merged_maf_blocks.block_ids[begin_pos] + 1,
                            0
                    );
                } else {
                    merged_block_id_intervals_tree.add(
                            merged_maf_blocks.block_ids[begin_pos],
                            merged_maf_blocks.block_ids[i - 1] + 1,
                            0
                    );
                }

                full_block_id_ranges += std::to_string(merged_maf_blocks.block_ids[begin_pos]);
                if ((i - 1) - begin_pos > 0) {
                    full_block_id_ranges += "-" + std::to_string(merged_maf_blocks.block_ids[i - 1]);
                }
                full_block_id_ranges += "_";

                begin_pos = i;
            }

            if (add_consensus) {
                is_block_in_a_merged_group[merged_maf_blocks.block_ids[i]] = true;
            }
        }

        if (inverted_merged_block_ids){
            merged_block_id_intervals_tree.add(
                    merged_maf_blocks.block_ids[merged_maf_blocks_size - 1],
                    merged_maf_blocks.block_ids[begin_pos] + 1,
                    0
            );
        } else {
            merged_block_id_intervals_tree.add(
                    merged_maf_blocks.block_ids[begin_pos],
                    merged_maf_blocks.block_ids[merged_maf_blocks_size - 1] + 1,
                    0
            );
        }
        full_block_id_ranges += std::to_string(merged_maf_blocks.block_ids[begin_pos]);
        if ((merged_maf_blocks_size - 1) - begin_pos > 0) {
            full_block_id_ranges += "-" + std::to_string(merged_maf_blocks.block_ids[merged_maf_blocks_size - 1]);
        }

        block_id_ranges_vector.push_back(block_id_ranges);
    }

    if (produce_maf) {
        bool contains_loops = false;

        auto maf = std::make_unique<ska::flat_hash_map<std::string, std::vector<maf_partial_row_t>>>();
        for (const auto &maf_prows : merged_maf_blocks.rows) {
            if (!contains_loops && maf_prows.second.size() > 1) {
                contains_loops = true;
            }

            maf->insert(std::pair<std::string, std::vector<maf_partial_row_t>>(
                    maf_prows.first,
                    std::vector<maf_partial_row_t>()
            ));

            for (auto &maf_prow : maf_prows.second) {
                (*maf)[maf_prows.first].push_back(
                        {
                                maf_prow.record_start,
                                maf_prow.seq_size,
                                maf_prow.is_reversed,
                                maf_prow.path_length,
                                maf_prow.aligned_seq
                        }
                );
            }
        }
        if (add_consensus) {
            uint64_t merged_consensus_seq_size = 0;
            uint64_t merged_consensus_path_length = 0;
            std::string merged_consensus_aligned_seq;

            uint64_t length_alignment = merged_maf_blocks.rows.begin()->second.begin()->aligned_seq.length();
            uint64_t start_cons_pos_in_alignment = 0;
            for (auto &maf_cons_prow : merged_maf_blocks.consensus_rows) {
                if (merged_maf_blocks_size == 1 || preserve_unmerged_consensus) {
                    std::string gapped_cons;

                    uint64_t pos;
                    for (pos = 0;
                         pos < start_cons_pos_in_alignment; pos++) { gapped_cons += "-"; }

                    gapped_cons += maf_cons_prow.second.aligned_seq;

                    pos += maf_cons_prow.second.aligned_seq.length();
                    for (; pos < length_alignment; pos++) { gapped_cons += "-"; }

                    maf->insert(std::pair<std::string, std::vector<maf_partial_row_t>>(
                            maf_cons_prow.first,
                            std::vector<maf_partial_row_t>()
                    ));

                    (*maf)[maf_cons_prow.first].push_back(
                            {
                                    maf_cons_prow.second.record_start,
                                    maf_cons_prow.second.seq_size,
                                    maf_cons_prow.second.is_reversed,
                                    maf_cons_prow.second.path_length,
                                    gapped_cons
                            }
                    );

                    clear_string(gapped_cons);

                    start_cons_pos_in_alignment += maf_cons_prow.second.aligned_seq.length();
                }

                // Manage merged consensus sequence
                if (merged_maf_blocks_size > 1) {
                    merged_consensus_seq_size += maf_cons_prow.second.seq_size;
                    merged_consensus_path_length += maf_cons_prow.second.path_length;
                    merged_consensus_aligned_seq += maf_cons_prow.second.aligned_seq;
                }
            }

            if (merged_maf_blocks_size > 1) {
                // Write the merged consensus

                maf->insert(std::pair<std::string, std::vector<maf_partial_row_t>>(
                        consensus_base_name + block_id_ranges + " ",
                        std::vector<maf_partial_row_t>()
                ));

                (*maf)[consensus_base_name + block_id_ranges + " "].push_back(
                        {
                                merged_maf_blocks.consensus_rows.begin()->second.record_start,
                                merged_consensus_seq_size,
                                merged_maf_blocks.consensus_rows.begin()->second.is_reversed,
                                merged_consensus_path_length,
                                merged_consensus_aligned_seq
                        }
                );
            }

            clear_string(merged_consensus_aligned_seq);
        }

        out_maf << "a blocks=" << full_block_id_ranges << " loops="
                << (contains_loops ? "true" : "false");
        if (merged_maf_blocks_size > 1) {
            out_maf << " merged=true";

            if (fraction_below_threshold) {
                out_maf << " below_thresh=true";
            }
        }
        out_maf << std::endl;

        write_maf_rows(out_maf, *maf);

    }

    // Cleaning
    /*
    clear_string(block_id_ranges);

    clear_vector(merged_maf_blocks.block_ids);
    for (const auto &maf_prows : merged_maf_blocks.rows) {
        for (auto &maf_prow : maf_prows.second) {
            clear_string(maf_prow.aligned_seq);
        }
    }
    merged_maf_blocks.rows.clear();

    for (auto &maf_cons_prow : merged_maf_blocks.consensus_rows) {
        clear_string(maf_cons_prow.first);
        clear_string(maf_cons_prow.second.aligned_seq);
    }
    clear_vector(merged_maf_blocks.consensus_rows);
    */
}

odgi::graph_t* smooth_and_lace(const xg::XG &graph,
                               blockset_t*& blockset,
                               int poa_m, int poa_n,
                               int poa_g, int poa_e,
                               int poa_q, int poa_c,
                               const bool& adaptive_poa_params,
                               const uint64_t &kmer_size,
                               float poa_padding_fraction,
                               uint64_t max_block_depth_for_padding_more,
                               bool local_alignment,
                               int n_threads,
                               int n_poa_threads,
                               std::string &path_output_maf, std::string &maf_header,
                               bool merge_blocks, bool preserve_unmerged_consensus, double contiguous_path_jaccard,
                               bool use_abpoa,
                               const std::string &consensus_base_name,
                               std::vector<std::string>& consensus_path_names,
                               bool write_fasta_blocks,
                               uint64_t max_merged_groups_in_memory) {

    bool add_consensus = !consensus_base_name.empty();

    //
    // record the start and end points of all the path ranges and the consensus
    //
    uint64_t block_count = blockset->size();
    auto _block_graphs = std::make_unique<std::vector<std::string*>>(block_count, nullptr);
    auto& block_graphs = *_block_graphs; // get a ref

    auto get_block_graph =
        [&](const uint64_t& block_id) {
            std::string data;
            zstdutil::DecompressString(*block_graphs[block_id], data);
            stringstream ss;
            ss << data;
            ss.seekg(0,std::ios_base::beg);
            auto block_graph = std::make_unique<odgi::graph_t>();
            block_graph->deserialize_members(ss);
            return block_graph;
        };

    auto save_block_graph =
        [&](const uint64_t& block_id,
            const odgi::graph_t* block_graph) {
            std::stringstream ss;
            block_graph->serialize_members(ss);
            std::string*& s = block_graphs[block_id];
            if (s == nullptr) {
                s = new std::string;
            } else {
                s->clear();
            }
            zstdutil::CompressString(ss.str(), *s);
        };

    // mapping from path fragments to block graphs
    auto _path_mapping_tmp = temp_file::create();
    auto path_mapping_ptr = std::make_unique<mmmulti::set<path_position_range_t>>(_path_mapping_tmp);
    auto& path_mapping = *path_mapping_ptr;
    path_mapping.open_writer();

    // mapping from block to consensus ids
    std::vector<path_handle_t> consensus_mapping(add_consensus ? blockset->size() : 0);

    std::vector<IITree<uint64_t, uint64_t>> merged_block_id_intervals_tree_vector;
    std::vector<std::string> block_id_ranges_vector;
    ska::flat_hash_set<uint64_t> inverted_merged_block_id_intervals_ranks; // IITree can't store inverted intervals

    std::atomic<uint64_t> num_flipped_graphs(0);
    atomicbitvector::atomic_bv_t blok_to_flip(blockset->size());

    std::vector<bool> is_block_in_a_merged_group((add_consensus && merge_blocks) ? blockset->size() : 0);

    {
        bool produce_maf = !path_output_maf.empty();

        // If merged consensus sequences have to be embedded, this structures are needed to keep the blocks' coordinates,
        // but the sequences will be considered (and kept in memory) only if a MAF has to be produced
        std::vector<std::unique_ptr<ska::flat_hash_map<std::string, std::vector<maf_partial_row_t>>>> mafs(produce_maf || (add_consensus && merge_blocks) ? blockset->size() : 0);
        atomicbitvector::atomic_bv_t mafs_ready(produce_maf || (add_consensus && merge_blocks) ? blockset->size() : 0);

        auto write_maf_lambda = [&]() {
            if (produce_maf || (add_consensus && merge_blocks)) {
                uint64_t block_id = 0;

                std::deque<std::unique_ptr<maf_t>> merged_maf_blocks_queue;

                uint64_t num_blocks = blockset->size();

                std::ofstream out_maf;

                if (produce_maf) {
                    out_maf.open(path_output_maf.c_str());
                    out_maf << maf_header << std::endl;
                }

                while (block_id < num_blocks) {
                    if (mafs_ready.test(block_id)) {
                        //std::cerr << "block_id (" << block_id << ")" << std::endl;

                        uint64_t num_seq_in_block = mafs[block_id]->size();
                        //for (auto const& path_to_maf_rows : mafs[block_id]) { num_seq_in_block += path_to_maf_rows.second.size(); }

                        bool merged = false;
                        bool prep_new_merge_group = false;
                        bool fraction_below_threshold = false;

                        std::string consensus_name;
                        if (add_consensus){
                            consensus_name = consensus_base_name + std::to_string(block_id);
                        }

                        bool flip_block_before_merging_in_the_group = false;
                        int64_t index_group_where_merge = -1;
                        int8_t new_block_on_the_left_in_the_group = - 1; // -1) undefined; 0) the new block go on the right; the new block go on the left

                        if (merge_blocks) {
                            if (merged_maf_blocks_queue.empty()) {
                                merged_maf_blocks_queue.push_back(std::make_unique<maf_t>());
                                index_group_where_merge = 0;

                                merged = true;
                            } else {
                                double best_jaccard = -1;

                                for(uint64_t num_group = 0; num_group < merged_maf_blocks_queue.size(); ++num_group) {
                                    auto& merged_maf_blocks = *merged_maf_blocks_queue[num_group];

                                    int8_t new_block_on_the_left = merged_maf_blocks.block_ids.size() > 1 ? (
                                            (merged_maf_blocks.block_ids.front() > merged_maf_blocks.block_ids.back() ? 1 : 0)
                                    ) : -1;

                                    bool flip_block_before_merging;

                                    uint64_t num_contiguous_ranges;
                                    for (auto flip_block : {false, true}) {
                                        flip_block_before_merging = flip_block;

                                        merged = true;
                                        num_contiguous_ranges = 0;

                                        for (auto const& path_to_maf_rows : *mafs[block_id]) {
                                            // Do not check the consensus (always forward)
                                            if (path_to_maf_rows.first != consensus_name){
                                                // To merge a block, it has to contains new sequences...
                                                if (merged_maf_blocks.rows.count(path_to_maf_rows.first) != 0) {
                                                    // ...or mergeable ones.

                                                    bool found_contiguous_row = false;
                                                    for (auto &maf_row : path_to_maf_rows.second) {
                                                        uint64_t maf_row_record_start = flip_block ? maf_row.path_length - (maf_row.record_start + maf_row.seq_size) : maf_row.record_start;

                                                        for (auto &merged_maf_prow : merged_maf_blocks.rows[path_to_maf_rows.first]) {
                                                            if ((flip_block ^ maf_row.is_reversed) == merged_maf_prow.is_reversed) {
                                                                if (flip_block ^ maf_row.is_reversed) {
                                                                    if ((merged_maf_prow.path_length - merged_maf_prow.record_start) ==
                                                                        (maf_row.path_length - (maf_row_record_start + maf_row.seq_size))) {
                                                                        // merged_maf_row_end == maf_row_begin, new row on the right

                                                                        if (new_block_on_the_left == -1 || new_block_on_the_left == 1) {
                                                                            // the row is the first one, or the merge continues on the left

                                                                            new_block_on_the_left = 1;

                                                                            found_contiguous_row = true;
                                                                            num_contiguous_ranges += 1;
                                                                            break;
                                                                        }
                                                                    } else if ((maf_row.path_length - maf_row_record_start) ==
                                                                               (merged_maf_prow.path_length - (merged_maf_prow.record_start + merged_maf_prow.seq_size))) {
                                                                        // maf_row_end == merged_maf_row_begin, new row on the left

                                                                        if (new_block_on_the_left == -1 || new_block_on_the_left == 0) {
                                                                            // the row is the first one, or the merge continues on the right

                                                                            new_block_on_the_left = 0;

                                                                            found_contiguous_row = true;
                                                                            num_contiguous_ranges += 1;
                                                                            break;
                                                                        }
                                                                    }
                                                                } else {
                                                                    if ((merged_maf_prow.record_start + merged_maf_prow.seq_size) == maf_row_record_start) {
                                                                        // merged_maf_row_end == maf_row_begin, new row on the right

                                                                        if (new_block_on_the_left == -1 || new_block_on_the_left == 0) {
                                                                            // the row is the first one, or the merge continues on the right

                                                                            new_block_on_the_left = 0;

                                                                            found_contiguous_row = true;
                                                                            num_contiguous_ranges += 1;
                                                                            break;
                                                                        }
                                                                    } else if ((maf_row_record_start + maf_row.seq_size) == merged_maf_prow.record_start) {
                                                                        // maf_row_end == merged_maf_row_begin, new row on the left

                                                                        if (new_block_on_the_left == -1 || new_block_on_the_left == 1) {
                                                                            // the row is the first one, or the merge continues on the left

                                                                            new_block_on_the_left = 1;

                                                                            found_contiguous_row = true;
                                                                            num_contiguous_ranges += 1;
                                                                            break;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }

                                                        // Commented out because we want to check if all ranges are mergeable
                                                        //if (found_contiguous_row) { break; }
                                                    }

                                                    if (!found_contiguous_row) {
                                                        merged = false; // Current block not mergeable
                                                        break;
                                                    }
                                                }
                                            }
                                        }

                                        if (merged) {
                                            uint64_t num_ranges_in_merged_block = 0;
                                            for (const auto &maf_prows : merged_maf_blocks.rows) { num_ranges_in_merged_block += maf_prows.second.size(); }

                                            uint64_t num_ranges_in_block_to_merge = 0;
                                            for (const auto &maf_prows : *mafs[block_id]) { num_ranges_in_block_to_merge += maf_prows.second.size(); }

                                            double current_contiguous_path_jaccard =
                                                    (double) num_contiguous_ranges /
                                                    (double) (num_ranges_in_block_to_merge - (add_consensus ? 1 : 0) + num_ranges_in_merged_block -
                                                            num_contiguous_ranges);

                                            if (current_contiguous_path_jaccard >= contiguous_path_jaccard && current_contiguous_path_jaccard > best_jaccard) {
                                                best_jaccard = current_contiguous_path_jaccard;

                                                flip_block_before_merging_in_the_group = flip_block_before_merging;
                                                index_group_where_merge = num_group;
                                                new_block_on_the_left_in_the_group = new_block_on_the_left;

                                                //// Greedy approach: both orientation are not always checked
                                                //break;
                                            }
                                        }
                                    }
                                }

                                fraction_below_threshold = best_jaccard > -1 && best_jaccard < contiguous_path_jaccard;
                            }

                            merged = index_group_where_merge > -1;
                        }
                        //std::cerr << "blockId " << block_id << " will be merged: " << (merged ? "yes" : "no") << std::endl;

                        // If mergeable...,
                        if (merged) {
                            // ..then merge

                            //std::cerr << "new_block_on_the_left_in_the_group " << (new_block_on_the_left_in_the_group == 1) << std::endl;

                            auto& merged_maf_blocks = merged_maf_blocks_queue[index_group_where_merge];

                            _put_block_in_group(
                                    *merged_maf_blocks, block_id, num_seq_in_block, consensus_name, mafs,
                                    new_block_on_the_left_in_the_group == 1,
                                    flip_block_before_merging_in_the_group);

                            //_clear_maf_block(mafs[block_id]);
                            mafs[block_id].reset(nullptr);

                            if (flip_block_before_merging_in_the_group) {
                                blok_to_flip.set(block_id);
                                ++num_flipped_graphs;
                            }
                        }

                        if (!merged) {
                            if (merged_maf_blocks_queue.size() >= max_merged_groups_in_memory){
                                // Write the merged group on the left
                                auto& merged_maf_blocks = merged_maf_blocks_queue.front();

                                _write_merged_maf_blocks(*merged_maf_blocks,
                                                         inverted_merged_block_id_intervals_ranks,
                                                         merged_block_id_intervals_tree_vector,
                                                         block_id_ranges_vector,
                                                         is_block_in_a_merged_group,
                                                         produce_maf, out_maf,
                                                         add_consensus, consensus_base_name,
                                                         fraction_below_threshold,
                                                         preserve_unmerged_consensus);

                                // Remove the merged group on the left from the queue
                                //delete merged_maf_blocks; //merged_maf_blocks_queue.front();
                                merged_maf_blocks_queue.pop_front();
                            }

                            // quietly groom by flipping the block to prefer the forward orientation of the lowest-ranked path
                            //

                            auto block_graph = get_block_graph(block_id);
                            uint64_t first_id = std::numeric_limits<uint64_t>::max();
                            path_handle_t groom_target_path;
                            block_graph->for_each_path_handle(
                                [&](const path_handle_t& p) {
                                    auto name_range = block_graph->get_path_name(p);
                                    auto path_name = name_range.substr(0, name_range.find_last_of('_'));
                                    uint64_t id = as_integer(graph.get_path_handle(path_name));
                                    if (id < first_id) {
                                        groom_target_path = p;
                                        first_id = id;
                                    }
                                });
                            assert(first_id < std::numeric_limits<uint64_t>::max());
                            bool flip_block = block_graph->get_is_reverse(
                                block_graph->get_handle_of_step(
                                    block_graph->path_begin(groom_target_path)));

                            // Put the current block in a new group on the right
                            merged_maf_blocks_queue.push_back(std::make_unique<maf_t>());
                            auto& merged_maf_blocks = merged_maf_blocks_queue.back();

                            // The last merge failed and the current un-merged block becomes the first one of the next group
                            _put_block_in_group(
                                    *merged_maf_blocks, block_id, num_seq_in_block, consensus_name, mafs,
                                    false,
                                    flip_block);

                            //_clear_maf_block(mafs[block_id]);
                            //delete mafs[block_id];
                            mafs[block_id].reset(nullptr);
                        }

                        block_id++;

                        /* TO REMOVE
                         if (!merged || is_last_block) {
                            if (!merged && !prep_new_merge_group) {
                                if (produce_maf) {
                                    bool contains_loops = false;
                                    ska::flat_hash_map<path_handle_t> seen_paths;
                                    for (auto &path_range : blockset->get_block(block_id).path_ranges) {
                                        path_handle_t path = graph.get_path_handle_of_step(path_range.begin);
                                        if (seen_paths.count(path)) {
                                            contains_loops = true;
                                            break;
                                        } else {
                                            seen_paths.insert(path);
                                        }
                                    }
                                    seen_paths.clear();

                                    out_maf << "a blocks=" + std::to_string(block_id) << " loops="
                                            << (contains_loops ? "true" : "false") << std::endl;
                                    write_maf_rows(out_maf, mafs[block_id]);
                                }

                                _clear_maf_block(mafs[block_id]);
                            }
                        }*/
                    }

                    std::this_thread::sleep_for(std::chrono::nanoseconds(1));
                }

                while (!merged_maf_blocks_queue.empty()) {
                    auto& merged_maf_blocks = merged_maf_blocks_queue.front();

                    _write_merged_maf_blocks(*merged_maf_blocks,
                                             inverted_merged_block_id_intervals_ranks,
                                             merged_block_id_intervals_tree_vector,
                                             block_id_ranges_vector,
                                             is_block_in_a_merged_group,
                                             produce_maf, out_maf,
                                             add_consensus, consensus_base_name,
                                             false,
                                             preserve_unmerged_consensus);

                    // Remove the merged group on the left from the queue
                    //delete merged_maf_blocks; //merged_maf_blocks_queue.front();
                    merged_maf_blocks_queue.pop_front();
                }

                if (produce_maf) {
                    out_maf.close();
                }

                //clear_vector(mafs);
            }
        };
        std::thread write_maf_thread(write_maf_lambda);

        std::stringstream poa_banner;
        poa_banner << "[smoothxg::smooth_and_lace] applying "
                   << (local_alignment ? "local" : "global") << " "
                   << (use_abpoa ? "abPOA" : "SPOA")
                   << " to " << blockset->size() << " blocks:";
        progress_meter::ProgressMeter poa_progress(blockset->size(), poa_banner.str());
        std::unique_ptr<ska::flat_hash_map<std::string, std::vector<maf_partial_row_t>>> empty_maf_block(nullptr);

        // Smooth blocks
#pragma omp parallel for schedule(dynamic,1) num_threads(n_poa_threads)
        for (uint64_t block_id = 0; block_id < blockset->size(); ++block_id) {
            auto block = blockset->get_block(block_id);

            std::string consensus_name;
            if (add_consensus){
                consensus_name = consensus_base_name + std::to_string(block_id);
            }

            // std::cerr << "on block " << block_id+1 << " of " << blockset->size() << std::endl;
            odgi::graph_t* block_graph = nullptr;
            if (produce_maf || (add_consensus && merge_blocks)) {
                mafs[block_id] = std::make_unique<ska::flat_hash_map<std::string, std::vector<maf_partial_row_t>>>();
            }

            uint64_t max_seq_len = 0;
            for (auto &path_range : block.path_ranges) {
                uint64_t seq_length = 0;
                for (step_handle_t step = path_range.begin; step != path_range.end;
                     step = graph.get_next_step(step)) {
                    const auto h = graph.get_handle_of_step(step);
                    const auto l = graph.get_length(h);
                    seq_length += l;
                }

                max_seq_len = std::max(max_seq_len, seq_length);
            }

            int poa_padding = 0;
            if (poa_padding_fraction > 0) {
                if (block.path_ranges.size() <= max_block_depth_for_padding_more) {
                    // min amount of flanking sequences to add
                    // avoids corner case of padding failure with shallow short blocks
                    poa_padding = 47;
                }
                // always pad a fraction of our sequence length
                // note that this can be very short in deep blocks that tend to have many short sequences
                // (this is because we establish blocks with a target mass and maximum length)
                {
                    float average_seq_len = 0.0;
                    // In blocks not too deep we can increase the padding size
                    for (auto &path_range : block.path_ranges) {
                        const path_handle_t path_handle = graph.get_path_handle_of_step(path_range.begin);
                        for (step_handle_t step = path_range.begin; step != path_range.end;
                        step = graph.get_next_step(step)) {
                            const auto h = graph.get_handle_of_step(step);
                            average_seq_len += (float)graph.get_length(h);
                        }
                    }
                    average_seq_len /= (float)block.path_ranges.size();
                    poa_padding = std::max((int)(average_seq_len * poa_padding_fraction), poa_padding);
                }
            }

            // Set/default penalties
            int poa_m_to_use = poa_m;
            int poa_n_to_use = poa_n;
            int poa_g_to_use = poa_g;
            int poa_e_to_use = poa_e;
            int poa_q_to_use = poa_q;
            int poa_c_to_use = poa_c;

            // Estimate the pairwise identity in the block for tuning the POA penalties
            // Avoid the identity estimation for too-deep blocks (todo random sampling for deep block???)
            if (adaptive_poa_params && block.path_ranges.size() > 1 && block.path_ranges.size() <= max_block_depth_for_padding_more) {
                // Deduplication (todo eventually, keep it simple for now)

                // Prepare sequences
                std::vector<std::string* > seqs; seqs.reserve(block.path_ranges.size());
                for (const auto& path_range : block.path_ranges) {
                    auto seq = new std::string();
                    for (step_handle_t step = path_range.begin; step != path_range.end; step = graph.get_next_step(step)) {
                        seq->append(graph.get_sequence(graph.get_handle_of_step(step)));
                    }

                    // We can't compute the hashes for what is shorter than the kmer size
                    // Skip too short sequences
                    if (seq->size() >= 8 * kmer_size) {
                        seqs.push_back(seq);
                    } else {
                        delete seq;
                    }
                }

                // Check if there are still sequences to compare
                if (seqs.size() > 1) {
                    // Compute hashes
                    std::vector<std::vector<mkmh::hash_t>> seq_hashes;
                    std::vector<int> seq_hash_lens;

                    seq_hashes.resize(seqs.size());
                    seq_hash_lens.resize(seqs.size());
                    rkmh::hash_sequences(seqs, seq_hashes, seq_hash_lens, kmer_size);

                    // All-vs-All comparison
                    std::vector<float> estimated_distances; //todo on-the-fly percentile computation to avoid a big vector in memory?

                    estimated_distances.reserve(seqs.size() * (seqs.size() - 1) / 2); // N * (N - 1) / 2 comparisons
                    for (uint64_t i = 0; i < seqs.size(); ++i) {
                        for (uint64_t j = i + 1; j < seqs.size(); ++j) {
                            const float est_identity = 1.0 - rkmh::compare(seq_hashes[i], seq_hashes[j], kmer_size, true);
                            estimated_distances.push_back(est_identity);
                        }
                    }

                    // Take 30% percentile as identity threshold (70% of the pairs have identity >= to this value)
                    std::sort(estimated_distances.begin(), estimated_distances.end());
                    const float est_identity_threshold = std::max((float)0.7, estimated_distances[(estimated_distances.size() - 1) * 0.30]);

                    // Tune POA penalties
                    if (est_identity_threshold >= 0.99) {
                        poa_m_to_use = 1;
                        poa_n_to_use = 19;
                        poa_g_to_use = 39;
                        poa_e_to_use = 3;
                        poa_q_to_use = 81;
                        poa_c_to_use = 1;
                    } else if (est_identity_threshold >= 0.98) {
                        poa_m_to_use = 1;
                        poa_n_to_use = 13;
                        poa_g_to_use = 31;
                        poa_e_to_use = 3;
                        poa_q_to_use = 51;
                        poa_c_to_use = 1;
                    } else if (est_identity_threshold >= 0.97) {
                        poa_m_to_use = 1;
                        poa_n_to_use = 9;
                        poa_g_to_use = 16;
                        poa_e_to_use = 2;
                        poa_q_to_use = 41;
                        poa_c_to_use = 1;
                    } else if (est_identity_threshold >= 0.95) {
                        poa_m_to_use = 1;
                        poa_n_to_use = 7;
                        poa_g_to_use = 11;
                        poa_e_to_use = 2;
                        poa_q_to_use = 33;
                        poa_c_to_use = 1;
                    } else if (est_identity_threshold >= 0.90) {
                        poa_m_to_use = 1;
                        poa_n_to_use = 4;
                        poa_g_to_use = 6;
                        poa_e_to_use = 2;
                        poa_q_to_use = 26;
                        poa_c_to_use = 1;
                    } // else use the set/default penalties
                }

                for (auto& seq : seqs) {
                    delete seq;
                }
            }

            if (use_abpoa) {
                block_graph = smooth_abpoa(graph,
                                           block,
                                           block_id,
                                           poa_m_to_use,
                                           poa_n_to_use,
                                           poa_g_to_use,
                                           poa_e_to_use,
                                           poa_q_to_use,
                                           poa_c_to_use,
                                           poa_padding,
                                           local_alignment,
                                           (produce_maf || (add_consensus && merge_blocks)) ? mafs[block_id] : empty_maf_block,
                                           produce_maf,
                                           true, // banded alignment
                                           consensus_name,
                                           write_fasta_blocks);
            } else {
                block_graph = smooth_spoa(graph,
                                          block,
                                          block_id,
                                          poa_m_to_use,
                                          -poa_n_to_use,
                                          -poa_g_to_use,
                                          -poa_e_to_use,
                                          -poa_q_to_use,
                                          -poa_c_to_use,
                                          poa_padding,
                                          local_alignment,
                                          (produce_maf || (add_consensus && merge_blocks)) ? mafs[block_id] : empty_maf_block,
                                          produce_maf,
                                          consensus_name,
                                          write_fasta_blocks);
            }

            // std::cerr << std::endl;
            // std::cerr << "After block graph. Exiting for now....." <<
            // std::endl; exit(0);
            if (block_graph->get_node_count() > 0) {
                // auto& block_graph = block_graphs.back();
                // record the start and end paths
                // nb: the path order is the same in the input block and output
                // graph
                uint64_t path_id = 0;
                for (auto &path_range : block.path_ranges) {
                    auto path_handle =
                        graph.get_path_handle_of_step(path_range.begin);
                    auto last_step = graph.get_previous_step(path_range.end);
                    path_mapping.append(
                        {path_handle, // target path
                         graph.get_position_of_step(
                             path_range.begin), // start position
                         (graph.get_position_of_step(
                             last_step) // end position
                          + graph.get_length(
                              graph.get_handle_of_step(last_step))),
                         as_path_handle(++path_id), block_id});
                }
                // make the graph

                // record the consensus path
                if (add_consensus) {
                    // record our consensus handle for later setup
                    consensus_mapping[block_id] = block_graph->get_path_handle(consensus_name);
                }
            }
            save_block_graph(block_id, block_graph);
            delete block_graph;
            poa_progress.increment(1);
            if (produce_maf || (add_consensus && merge_blocks)){
                mafs_ready.set(block_id);
            }
        }

        poa_progress.finish();

        write_maf_thread.join();
    }

    // Flip graphs
    {
        std::stringstream flip_graphs_banner;
        flip_graphs_banner << "[smoothxg::smooth_and_lace] flipping " << num_flipped_graphs << " block graphs:";
        progress_meter::ProgressMeter flip_graphs_progress(block_count, flip_graphs_banner.str());

#pragma omp parallel for schedule(dynamic,1)
        for (uint64_t block_id = 0; block_id < block_count; ++block_id) {
            if (blok_to_flip.test(block_id)) {
                auto* flipped_graph = new odgi::graph_t();

                ska::flat_hash_map<handle_t, handle_t> forward_translation;

                auto block_graph = get_block_graph(block_id);

                // make the flipped nodes
                block_graph->for_each_handle(
                        [&](const handle_t& old_handle) {
                            handle_t rev_handle = flipped_graph->create_handle(
                                    block_graph->get_sequence(block_graph->flip(old_handle)),
                                    block_graph->get_id(old_handle));
                            forward_translation[old_handle] = rev_handle;
                        });

                // make the flipped edges
                block_graph->for_each_edge([&](const edge_t& edge) {
                    // get the two sides in the correct orientation
                    handle_t rev_left, rev_right;
                    if (block_graph->get_is_reverse(edge.first)) {
                        rev_left = flipped_graph->flip(forward_translation[block_graph->flip(edge.first)]);
                    }
                    else {
                        rev_left = flipped_graph->flip(forward_translation[edge.first]);
                    }

                    if (!block_graph->get_is_reverse(edge.second)) {
                        rev_right = flipped_graph->flip(forward_translation[block_graph->flip(edge.second)]);
                    }
                    else {
                        rev_right = flipped_graph->flip(forward_translation[edge.second]);
                    }

                    // actually make the edge
                    flipped_graph->create_edge(rev_left, rev_right);
                });

                std::string consensus_name;
                if (add_consensus){
                    consensus_name = consensus_base_name + std::to_string(block_id);
                }

                block_graph->for_each_path_handle(
                        [&](const path_handle_t& old_path) {
                            string path_name = block_graph->get_path_name(old_path);
                            path_handle_t new_path = flipped_graph->create_path_handle(path_name);
                            if (path_name != consensus_name) {
                                // flip the path, but preserving its sequence
                                block_graph->for_each_step_in_path(old_path, [&](const step_handle_t& step) {
                                    handle_t old_handle = block_graph->get_handle_of_step(step);
                                    handle_t new_handle = flipped_graph->get_handle(
                                            block_graph->get_id(old_handle),
                                            !block_graph->get_is_reverse(old_handle));
                                    flipped_graph->append_step(new_path, new_handle);
                                });
                            } else {
                                // the consensus has to be encoded in reverse complement, but it remains in the forward strand
                                block_graph->for_each_step_in_path(old_path, [&](const step_handle_t& step) {
                                    handle_t old_handle = block_graph->get_handle_of_step(step);
                                    handle_t new_handle = flipped_graph->get_handle(
                                            block_graph->get_id(old_handle),
                                            block_graph->get_is_reverse(old_handle));
                                    flipped_graph->prepend_step(new_path, new_handle);
                                });
                            }

                        });

                save_block_graph(block_id, flipped_graph);
                delete flipped_graph;

                flip_graphs_progress.increment(1);
            }
        }
        flip_graphs_progress.finish();
    }

    std::cerr << "[smoothxg::smooth_and_lace] sorting path_mappings" << std::endl;
    // sort the path range mappings by path handle id, then start position
    // this will allow us to walk through them in order
    /*
    ips4o::parallel::sort(
        path_mapping.begin(), path_mapping.end(),
        [](const path_position_range_t &a, const path_position_range_t &b) {
            auto &a_id = as_integer(get_base_path(a));
            auto &b_id = as_integer(get_base_path(b));
            return (a_id < b_id || a_id == b_id && get_start_pos(a) < get_start_pos(b));
        });
    */
    path_mapping.index(n_threads);

    // build the sequence and edges into the output graph
    auto* smoothed = new odgi::graph_t();
    std::vector<path_handle_t> paths; // for parallel iteration

    // add the nodes and edges to the graph
    std::vector<uint64_t> id_mapping;
    {
        std::vector<std::unique_ptr<odgi::graph_t>> graphs(block_count);
#pragma omp parallel for schedule(dynamic,1)
        for (uint64_t idx = 0; idx < block_count; ++idx) {
            std::string data;
            zstdutil::DecompressString(*block_graphs[idx], data);
            stringstream ss;
            ss << data;
            ss.seekg(0,std::ios_base::beg);
            graphs[idx] = std::make_unique<odgi::graph_t>();
            graphs[idx]->deserialize_members(ss);
            delete block_graphs[idx];
        }
        _block_graphs.reset(nullptr); // we've decompressed these, now clear our block graphs

        std::stringstream add_graph_banner;
        add_graph_banner << "[smoothxg::smooth_and_lace] adding nodes from " << block_count << " graphs:";
        progress_meter::ProgressMeter add_graph_progress(block_count, add_graph_banner.str());

        for (uint64_t idx = 0; idx < block_count; ++idx) {
            uint64_t id_trans = smoothed->get_node_count();
            // record the id translation
            auto& block = graphs[idx];
            id_mapping.push_back(id_trans);
            if (block->get_node_count() == 0) {
                continue;
            }
            block->for_each_handle([&](const handle_t &h) {
                smoothed->create_handle(block->get_sequence(h));
            });
            add_graph_progress.increment(1);
        }

        std::stringstream add_edges_banner;
        add_edges_banner << "[smoothxg::smooth_and_lace] adding edges from " << block_count << " graphs:";
        progress_meter::ProgressMeter add_edges_progress(block_count, add_edges_banner.str());
        for (uint64_t idx = 0; idx < block_count; ++idx) {
            auto& id_trans = id_mapping[idx];
            auto& block = graphs[idx];
            block->for_each_edge([&](const edge_t &e) {
                smoothed->create_edge(
                        smoothed->get_handle(id_trans + block->get_id(e.first)),
                        smoothed->get_handle(id_trans + block->get_id(e.second)));
            });
            add_edges_progress.increment(1);
        }
        add_edges_progress.finish();

        // then for each path, ensure that it's embedded in the graph by walking through
        // its block segments in order and linking them up in the output graph
        std::stringstream lace_banner;
        lace_banner << "[smoothxg::smooth_and_lace] embedding " << path_mapping.size() << " path fragments:";
        progress_meter::ProgressMeter lace_progress(path_mapping.size(), lace_banner.str());
        for (uint64_t i = 0; i < path_mapping.size(); ++i) {
            path_position_range_t pos_range = path_mapping.read_value(i);
            step_handle_t last_step = {0, 0};
            bool first = true;
            uint64_t last_end_pos = 0;
            // add the path to the graph

            path_handle_t smoothed_path = smoothed->create_path_handle(
                graph.get_path_name(get_base_path(pos_range)));
            // walk the path from start to end
            while (true) {
                // if we find a segment that's not included in any block, we'll add
                // it to the final graph and link it in to do so, we detect a gap in
                // length, collect the sequence in the gap and add it to the graph
                // as a node then add it as a traversal to the path
                if (get_start_pos(pos_range) - last_end_pos > 0) {
                    assert(false); // assert that we've included all sequence in blocks
                }
                // write the path steps into the graph using the id translation
                auto block_id = get_block_id(pos_range);
                auto& block = graphs[block_id];
                auto id_trans = id_mapping.at(block_id);
                block->for_each_step_in_path(
                    get_target_path(pos_range), [&](const step_handle_t &step) {
                        handle_t h = block->get_handle_of_step(step);
                        handle_t t = smoothed->get_handle(block->get_id(h) + id_trans,
                                                          block->get_is_reverse(h));
                        smoothed->append_step(smoothed_path, t);
                        if (first) {
                            first = false;
                            // create edge between last and curr
                            if (as_integers(last_step)[0] != 0) {
                                smoothed->create_edge(
                                    smoothed->get_handle_of_step(last_step), t);
                            }
                        }
                    });
                last_step = smoothed->path_back(smoothed_path);
                last_end_pos = get_end_pos(pos_range);
                if (i + 1 == path_mapping.size() ||
                    get_base_path(path_mapping.read_value(i + 1)) != get_base_path(pos_range)) {
                    break;
                } else {
                    ++i;
                    pos_range = path_mapping.read_value(i);
                }
                lace_progress.increment(1);
            }
            // now add in any final sequence in the path
            // and add it to the path, add the edge
            if (graph.get_path_length(get_base_path(pos_range)) > last_end_pos) {
                assert(false); // assert that we've included all sequence in the blocks
            }
        }
        lace_progress.finish();

        path_mapping.close_reader();
        std::remove(_path_mapping_tmp.c_str());
        path_mapping_ptr.reset(nullptr);

        // now verify that smoothed has paths that are equal to the base graph
        // and that all the paths are fully embedded in the graph
        smoothed->for_each_path_handle(
            [&](const path_handle_t &path) {
                paths.push_back(path);
            });

        {
            std::stringstream validate_banner;
            validate_banner << "[smoothxg::smooth_and_lace] validating " << paths.size() << " path sequences:";
            progress_meter::ProgressMeter validate_progress(paths.size(), validate_banner.str());

#pragma omp parallel for schedule(dynamic,1)
            for (uint64_t i = 0; i < paths.size(); ++i) {
                auto path = paths[i];

                std::string orig_seq, smoothed_seq;
                graph.for_each_step_in_path(
                    graph.get_path_handle(smoothed->get_path_name(path)),
                    [&](const step_handle_t &step) {
                        orig_seq.append(graph.get_sequence(graph.get_handle_of_step(step)));
                    });
                smoothed->for_each_step_in_path(
                    path,
                    [&](const step_handle_t &step) {
                        smoothed_seq.append(smoothed->get_sequence(smoothed->get_handle_of_step(step)));
                    });
                if (orig_seq != smoothed_seq) {
                    std::cerr << "[smoothxg] error! path "
                              << smoothed->get_path_name(path)
                              << " was corrupted in the smoothed graph" << std::endl
                              << "original\t" << orig_seq << std::endl
                              << "smoothed\t" << smoothed_seq << std::endl;
                    exit(1);
                }

                validate_progress.increment(1);
            }
            validate_progress.finish();
        }

        if (!consensus_mapping.empty()) {
            std::cerr << "[smoothxg::smooth_and_lace] sorting consensus" << std::endl;

            // consensus path and connections

            // by definition, the consensus paths are embedded in our blocks, which simplifies
            // things we'll still need to add a new path for each consensus path

            // flag the blocks that we should include unmerged
            atomicbitvector::atomic_bv_t exclude_unmerged_consensus(block_count);

            // Is there something merged?
            if (!merged_block_id_intervals_tree_vector.empty()) {
#pragma omp parallel for schedule(dynamic,1)
                for (auto& merged_block_id_intervals_tree : merged_block_id_intervals_tree_vector) {
                    merged_block_id_intervals_tree.index();
                }

                if (!preserve_unmerged_consensus) {
                    std::cerr << "[smoothxg::smooth_and_lace] embedding consensus: removing redundant single consensus" << std::endl;

#pragma omp parallel for schedule(dynamic,1)
                    for (uint64_t id = 0; id < consensus_mapping.size(); ++id) {
                        if (is_block_in_a_merged_group[id]) {
                            exclude_unmerged_consensus.set(id);
                        }
                    }
                }
            }

            // all raw consensus paths
            std::vector<path_handle_t> consensus_paths(block_count);

            // Unmerged consensus sequences
            // First, create the path handles
            std::cerr << "[smoothxg::smooth_and_lace] embedding consensus: creating path handles" << std::endl;
            for (uint64_t id = 0; id < consensus_mapping.size(); ++id) {
                //for (auto &pos_range : consensus_mapping) {
                if (!exclude_unmerged_consensus.test(id)) {
                    auto& block = graphs[id];
                    consensus_paths[id] = smoothed->create_path_handle(
                            block->get_path_name(consensus_mapping[id]));
                } // else skip the embedding of the single consensus sequences
            }

            // Next, add the steps
            std::cerr << "[smoothxg::smooth_and_lace] embedding consensus: creating step handles" << std::endl;
#pragma omp parallel for schedule(dynamic,1)
            for (uint64_t id = 0; id < consensus_mapping.size(); ++id) {
                //for(auto& pos_range : consensus_mapping){
                if (exclude_unmerged_consensus.test(id)) {
                    continue; // skip the embedding for the single consensus sequence
                }
                auto& block = graphs[id];
                path_handle_t smoothed_path = consensus_paths[id];
                auto &id_trans = id_mapping[id];
                block->for_each_step_in_path(consensus_mapping[id], [&](const step_handle_t &step) {
                        handle_t h = block->get_handle_of_step(step);
                        handle_t t = smoothed->get_handle(block->get_id(h) + id_trans, block->get_is_reverse(h));
                        smoothed->append_step(smoothed_path, t);
                        // nb: by definition of our construction of smoothed
                        // the consensus paths should have all their edges embedded
                    });
            }

            // Merged consensus sequences
            if (!merged_block_id_intervals_tree_vector.empty()) {
                // First, create the path handles
                std::cerr << "[smoothxg::smooth_and_lace] embedding merged consensus: creating path handles" << std::endl;
                std::vector<path_handle_t> merged_consensus_paths;

                for (auto &block_id_ranges : block_id_ranges_vector) {
                    assert(!smoothed->has_path(consensus_base_name + block_id_ranges));
                    merged_consensus_paths.push_back(
                        smoothed->create_path_handle(consensus_base_name + block_id_ranges)
                        );
                }

                // Next, add the steps
                std::cerr << "[smoothxg::smooth_and_lace] embedding merged consensus: creating step handles" << std::endl;
                std::mutex consensus_path_is_merged_mutex;
                ska::flat_hash_set<uint64_t> consensus_path_is_merged;
                assert(merged_block_id_intervals_tree_vector.size() == block_id_ranges_vector.size());

#pragma omp parallel for schedule(dynamic,1)
                for (uint64_t i = 0; i < merged_block_id_intervals_tree_vector.size(); ++i) {
                    auto &merged_block_id_intervals_tree = merged_block_id_intervals_tree_vector[i];

                    bool inverted_intervals = inverted_merged_block_id_intervals_ranks.count(i) != 0;
                    path_handle_t consensus_path = merged_consensus_paths[i];

                    std::vector<size_t> merged_block_id_intervals;
                    merged_block_id_intervals_tree.overlap(0, block_count, merged_block_id_intervals);

                    uint64_t start_interval = 0;
                    uint64_t end_interval = merged_block_id_intervals.size() - 1;
                    int8_t step_interval = 1;
                    int8_t step = 1;
                    if (inverted_intervals) {
                        start_interval = merged_block_id_intervals.size() - 1;
                        end_interval = 0;
                        step_interval = -1;
                        step = -1;
                    }

                    for (uint64_t j = start_interval; j != (end_interval + step_interval); j += step_interval) {
                        auto &merged_block_id_interval_idx = merged_block_id_intervals[j];

                        uint64_t start = merged_block_id_intervals_tree.start(merged_block_id_interval_idx);
                        uint64_t end = merged_block_id_intervals_tree.end(merged_block_id_interval_idx) - 1;
                        if (inverted_intervals){
                            uint64_t tmp = start;
                            start = end;
                            end = tmp;

                            /*{
                              std::lock_guard<std::mutex> guard(consensus_path_is_merged_mutex);

                              std::cerr << i << ": start-end " << start << "-" << end <<std::endl;
                              }*/
                        }

                        for (uint64_t block_id = start; block_id != (end + step); block_id += step) {
                            {
                                std::lock_guard<std::mutex> guard(consensus_path_is_merged_mutex);

                                consensus_path_is_merged.insert(as_integer(consensus_paths[block_id]));
                            }

                            auto& block = graphs[block_id];
                            auto& id_trans = id_mapping[block_id];
                            block->for_each_step_in_path(
                                consensus_mapping[block_id],
                                [&](const step_handle_t &step) {
                                    handle_t h = block->get_handle_of_step(step);
                                    handle_t t = smoothed->get_handle(block->get_id(h) + id_trans, block->get_is_reverse(h));
                                    smoothed->append_step(consensus_path, t);
                                });
                        }
                    }

                    clear_string(block_id_ranges_vector[i]);
                }

                // now for each consensus path that's not been merged, and for each merged consensus path...
                // record our path handles for later use in consensus graph generation

                consensus_paths.erase(
                    std::remove_if(
                        consensus_paths.begin(), consensus_paths.end(),
                        [&consensus_path_is_merged](const path_handle_t& path) {
                            return consensus_path_is_merged.count(as_integer(path)) > 0;
                        }),
                    consensus_paths.end());

                consensus_paths.reserve(
                    consensus_paths.size()
                    + std::distance(merged_consensus_paths.begin(),
                                    merged_consensus_paths.end()));
                consensus_paths.insert(
                    consensus_paths.end(),
                    merged_consensus_paths.begin(),
                    merged_consensus_paths.end());

            }

            // todo: validate the consensus paths as well

            consensus_path_names.reserve(consensus_paths.size());
            for (auto &path : consensus_paths) {
                consensus_path_names.push_back(smoothed->get_path_name(path));
            }
        }
    }

    {
        std::stringstream embed_banner;
        embed_banner << "[smoothxg::smooth_and_lace] walking edges in "
                     << paths.size() << " paths:";
        progress_meter::ProgressMeter embed_progress(paths.size(), embed_banner.str());
        // embed all paths in the graph to ensure validity
        smoothed->for_each_path_handle(
            [&](const path_handle_t& path) {
                handle_t last;
                step_handle_t begin_step = smoothed->path_begin(path);
                smoothed->for_each_step_in_path(
                    path,
                    [&](const step_handle_t &step) {
                        handle_t h = smoothed->get_handle_of_step(step);
                        if (step != begin_step) {
                            smoothed->create_edge(last, h);
                        }
                        last = h;
                    });
                embed_progress.increment(1);
            });
        embed_progress.finish();
    }

    return smoothed;
}

void build_odgi_abPOA(abpoa_t *ab, abpoa_para_t *abpt, odgi::graph_t* output,
                      const std::vector<std::string> &sequence_names,
                      const std::vector<bool>& aln_is_reverse,
                      const std::string &consensus_name,
                      const int &padding_len,
                      bool include_consensus) {
    abpoa_seq_t *abs = ab->abs; abpoa_graph_t *abg = ab->abg;
    // how would this happen, and can we manage the error externally?
    if (abg->node_n <= 2) return;

    // traverse graph
    int *in_degree = (int*)_err_malloc(abg->node_n * sizeof(int));
    int n_seq = abs->n_seq;
    int **read_paths = (int**)_err_malloc(n_seq * sizeof(int*)), *read_path_i = (int*)_err_calloc(n_seq, sizeof(int));
    int i, j, cur_id, pre_id, out_id, *id;
    for (i = 0; i < abg->node_n; ++i)
        in_degree[i] = abg->node[i].in_edge_n;
    for (i = 0; i < n_seq; ++i)
        read_paths[i] = (int*)_err_malloc(abg->node_n * sizeof(int));

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
                //fprintf(stdout, "S\t%d\t%c\n", cur_id-1, ab_char256_table[abg->node[cur_id].base]);
                // add node to output graph
                std::string seq = std::string(1, ab_char256_table[abg->node[cur_id].base]);
                // std::cerr << "seq: " << seq << std::endl;
                output->create_handle(seq, cur_id-1);
                // output all links based pre_ids
                for (i = 0; i < abg->node[cur_id].in_edge_n; ++i) {
                    pre_id = abg->node[cur_id].in_id[i];
                    if (pre_id != ABPOA_SRC_NODE_ID){
                        //fprintf(stdout, "L\t%d\t+\t%d\t+\t0M\n", pre_id-1, cur_id-1);
                        output->create_edge(output->get_handle(pre_id-1), output->get_handle(cur_id-1));
                    }
                }
                // add node id to read path
                int b, read_id; uint64_t num, tmp;
                b = 0;
                for (i = 0; i < abg->node[cur_id].read_ids_n; ++i) {
                    for (j = 0; j < abg->node[cur_id].out_edge_n; ++j) {
                        num = abg->node[cur_id].read_ids[j][i];
                        while (num) {
                            tmp = num & -num;
                            read_id = ilog2_64(tmp);
                            read_paths[b+read_id][read_path_i[b+read_id]++] = cur_id-1;
                            num ^= tmp;
                        }
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
    for (i = 0; i < n_seq; ++i) {
        //fprintf(stdout, "P\t%s\t", sequence_names[i]);
        path_handle_t p = output->create_path_handle(sequence_names[i]);

        if (aln_is_reverse[i]) {
            for (j = read_path_i[i] - 1 - padding_len; j >= padding_len; --j) {
                //fprintf(stdout, "%d-", read_paths[i][j]);
                //if (j != 0) fprintf(stdout, ","); else fprintf(stdout, "\t*\n");
                output->append_step(p, output->flip(output->get_handle(read_paths[i][j])));
            }
        } else {
            for (j = padding_len; j < read_path_i[i] - padding_len; ++j) {
                //fprintf(stdout, "%d+", read_paths[i][j]);
                //if (j != read_path_i[i]-1) fprintf(stdout, ","); else fprintf(stdout, "\t*\n");
                output->append_step(p, output->get_handle(read_paths[i][j]));
            }
        }
    }
    if (include_consensus) {
        abpoa_cons_t *abc = ab->abc;
        int cons_i = 0; // Only the first consensus

        path_handle_t p = output->create_path_handle(consensus_name);

        //fprintf(stdout, "P\tConsensus_sequence"); if (abc->n_cons > 1) fprintf(stdout, "_%d", cons_i+1); fprintf(out_fp, "\t");
        for (i = 1; i < abc->cons_len[cons_i]; ++i) {
            cur_id = abc->cons_node_ids[cons_i][i];
            const uint64_t step_count = output->steps_of_handle(output->get_handle(cur_id-1)).size();
            if (step_count > 0) {
                // It is an handle supported by at least one original path too
                output->append_step(p, output->get_handle(cur_id-1));
            }
            //fprintf(stdout, "%d+", cur_id-1); if (i != abc->cons_len[cons_i]-1) fprintf(stdout, ","); else fprintf(stdout, "\t*\n");
        }
    }

    free(in_degree);
    for (i = 0; i < n_seq; ++i)
        free(read_paths[i]);
    free(read_paths);
    free(read_path_i);

    // Remove unused edges
    {
        std::vector<edge_t> edges_to_drop = odgi::algorithms::find_edges_exceeding_depth_limits(*output, 1, std::numeric_limits<uint64_t>::max());
        for (auto& edge : edges_to_drop) {
            output->destroy_edge(edge);
        }
    }

    // Remove unused nodes
    {
        std::vector<handle_t> handles_to_drop = odgi::algorithms::find_handles_exceeding_depth_limits(*output, 1, std::numeric_limits<uint64_t>::max());
        for (auto& handle : handles_to_drop) {
            output->destroy_handle(handle);
        }
    }
}

void build_odgi_SPOA(spoa::Graph& graph, odgi::graph_t* output,
                     const std::vector<std::string> &sequence_names,
                     const int &padding_len,
                     const std::vector<bool> &aln_is_reverse,
                     const std::string &consensus_name, bool include_consensus) {

    auto &nodes = graph.nodes();

    for (const auto& it : graph.nodes()) {
        std::string s(1, static_cast<char>(graph.decoder(it->code)));
        auto h = output->create_handle(s, it->id + 1);
    }
    for (const auto& it : graph.nodes()) {
        auto h = output->get_handle(it->id + 1);
        for (const auto& jt : it->outedges) {
            auto t = output->get_handle(jt->head->id + 1);
            output->create_edge(h, t);
        }
    }

    for (std::uint32_t i = 0; i < graph.sequences().size(); ++i) {
        path_handle_t p = output->create_path_handle(sequence_names[i]);
        std::vector<handle_t> steps;
        auto curr = graph.sequences()[i];
        while (true) {
            steps.emplace_back(output->get_handle(curr->id + 1));
            if (!(curr = curr->Successor(i))) {
                break;
            }
        }
        steps = std::vector<handle_t>(steps.begin() + padding_len, steps.end() - padding_len);
        if (aln_is_reverse[i]) {
            for (auto handle_itr = steps.rbegin(); handle_itr != steps.rend(); ++handle_itr) {
                output->append_step(p, output->flip(*handle_itr));
            }
        } else {
            for (auto &handle : steps) {
                output->append_step(p, handle);
            }
        }
    }

    if (include_consensus) {
        path_handle_t p = output->create_path_handle(consensus_name);
        for (std::uint32_t i = 0; i < graph.consensus().size(); ++i) {
            output->append_step(p, output->get_handle(graph.consensus()[i]->id + 1));
            /*
              // not clear if this is still needed
            const uint64_t step_count = output->steps_of_handle(output->get_handle(id + 1)).size();
            if (step_count > 0) {
                // It is an handle supported by at least one original path too
                output->append_step(p, output->get_handle(id + 1));
            }
            */
        }
    }

    // Remove unused edges
    {
        std::vector<edge_t> edges_to_drop = odgi::algorithms::find_edges_exceeding_depth_limits(*output, 1, std::numeric_limits<uint64_t>::max());
        for (auto& edge : edges_to_drop) {
            output->destroy_edge(edge);
        }
    }

    // Remove unused nodes
    {
        std::vector<handle_t> handles_to_drop = odgi::algorithms::find_handles_exceeding_depth_limits(*output, 1, std::numeric_limits<uint64_t>::max());
        for (auto& handle : handles_to_drop) {
            output->destroy_handle(handle);
        }
    }
}

} // namespace smoothxg
