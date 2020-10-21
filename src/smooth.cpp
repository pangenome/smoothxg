#include "smooth.hpp"
#include <cstring>

#include "deps/abPOA/src/seq.h"
#include "deps/abPOA/src/abpoa_graph.h"

#include "maf.hpp"

#include "progress.hpp"

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

void _clear_maf_block(uint64_t block_id, std::vector<std::vector<maf_row_t>> &mafs){
    for (auto& maf_row : mafs[block_id]){
        clear_string(maf_row.path_name);
        clear_string(maf_row.aligned_seq);
    }
    clear_vector(mafs[block_id]);
}

odgi::graph_t smooth_abpoa(const xg::XG &graph, const block_t &block, const uint64_t &block_id,
                           int poa_m, int poa_n, int poa_g,
                           int poa_e, int poa_q, int poa_c,
                           bool local_alignment,
                           std::vector<maf_row_t> *maf, bool keep_sequence,
                           bool banded_alignment,
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
    if (!banded_alignment) abpt->wb = -1;
    //abpt->zdrop = 100; // could be useful in local mode
    //abpt->end_bonus = 100; // also useful in local mode
    abpt->rev_cigar = 0;
    abpt->out_gfa = 1; // must be set to get the graph
    abpt->out_msa = maf != nullptr ? 1 : 0; // must be set when we extract the MSA
    abpt->out_cons = generate_consensus;
    abpt->amb_strand = 1; // we'll align both ways and check which is better
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
                nst_nt4_table[(int)seqs_[i][j]]; // TODO we make a c_str for every char in the string
        }
    }

    // variables to store result
    uint8_t **cons_seq; int **cons_cov, *cons_l, cons_n=0;
    uint8_t **msa_seq; int msa_l=0;

    int i, tot_n = n_seqs;
    uint8_t *is_rc = (uint8_t *)_err_malloc(n_seqs * sizeof(uint8_t));
    abpoa_reset_graph(ab, abpt, seq_lens[0]);

    std::vector<bool> aln_is_reverse;
    for (i = 0; i < n_seqs; ++i) {
        abpoa_res_t res;
        res.graph_cigar = 0, res.n_cigar = 0, res.is_rc = 0;
        res.traceback_ok = 1;
        abpt->rev_cigar = 0;
        bool aligned = -1 != abpoa_align_sequence_to_graph(ab, abpt, bseqs[i], seq_lens[i], &res);
        // nb: we should check if we should do anything special when !res->traceback_ok
        abpoa_add_graph_alignment(ab, abpt, bseqs[i], seq_lens[i], res, i, n_seqs);
        is_rc[i] = res.is_rc;
        if (res.is_rc) {
            aln_is_reverse.push_back(true);
        } else {
            aln_is_reverse.push_back(false);
        }
        if (res.n_cigar) {
            free(res.graph_cigar);
        }
    }

    if (maf != nullptr){
        abpoa_generate_rc_msa(ab, abpt, NULL, is_rc, tot_n, NULL, &msa_seq, &msa_l);
    }

   /*fprintf(stdout, ">Multiple_sequence_alignment\n");
   for (i = 0; i < n_seqs; ++i) {
       for (int j = 0; j < msa_l; ++j) {
           fprintf(stdout, "%c", "ACGTN-"[msa_seq[i][j]]);
       }
       fprintf(stdout, "\n");
   }*/

    if (generate_consensus) {
        abpoa_generate_consensus(ab, abpt, tot_n, NULL, &cons_seq, &cons_cov, &cons_l, &cons_n);
        if (ab->abg->is_called_cons == 0) {
            // TODO Abort mission here?
            err_printf("ERROR: no consensus sequence generated.\n");
            exit(1);
        }
        aln_is_reverse.push_back(false);

        /*fprintf(stdout, "=== output to variables ===\n");
        for (int i = 0; i < cons_n; ++i) {
            fprintf(stdout, ">Consensus_sequence\n");
            for (int j = 0; j < cons_l[i]; ++j)
                fprintf(stdout, "%c", nst_nt256_table[cons_seq[i][j]]);
            fprintf(stdout, "\n");
        }*/
    }

    if (maf != nullptr) {
        uint64_t num_seqs = n_seqs + (generate_consensus ? 1 : 0);
        for(uint64_t seq_rank = 0; seq_rank < num_seqs; seq_rank++) {
            std::basic_string<char> aligned_seq;

            std::string path_name;
            uint64_t seq_size;
            uint64_t path_length;
            uint64_t record_start;
            if (!generate_consensus || seq_rank < num_seqs - 1) {
                if (keep_sequence){
                    for (int j = 0; j < msa_l; ++j) {
                        aligned_seq += "ACGTN-"[msa_seq[seq_rank][j]];
                    }
                }

                auto path_handle = graph.get_path_handle_of_step(block.path_ranges[seq_rank].begin);

                path_name = graph.get_path_name(path_handle);
                path_length = graph.get_path_length(path_handle);

                // If the strand field is "-" then this is the start relative to the reverse-complemented source sequence
                uint64_t path_range_begin = graph.get_position_of_step(block.path_ranges[seq_rank].begin);
                auto last_step = graph.get_previous_step(block.path_ranges[seq_rank].end);
                record_start = aln_is_reverse[seq_rank] ?
                               (path_length - graph.get_position_of_step(last_step) -
                                graph.get_length(graph.get_handle_of_step(last_step))) :
                               path_range_begin;

                seq_size = seqs[seq_rank].size(); // <==> block.path_ranges[seq_rank].length
            } else {
                // The last sequence is the gapped consensus

                if (keep_sequence){
                    int j, k, aligned_id, rank;
                    i = ab->abg->node[ABPOA_SRC_NODE_ID].max_out_id;
                    int last_rank = 1;
                    while (i != ABPOA_SINK_NODE_ID) {
                        rank = abpoa_graph_node_id_to_msa_rank(ab->abg, i);
                        for (k = 0; k < ab->abg->node[i].aligned_node_n; ++k) {
                            aligned_id = ab->abg->node[i].aligned_node_id[k];
                            rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(ab->abg, aligned_id));
                        }
                        // last_rank -> rank : -
                        for (k = last_rank; k < rank; ++k) aligned_seq += '-';
                        // rank : base
                        aligned_seq += "ACGTN"[ab->abg->node[i].base];
                        last_rank = rank+1;
                        i = ab->abg->node[i].max_out_id;
                    }
                    // last_rank -> msa_l:-
                    for (k = last_rank; k <= msa_l; ++k) aligned_seq += '-';
                }

                path_name = consensus_name;
                path_length = cons_l[0];
                record_start = 0;
                seq_size = path_length;
            }

            maf->push_back({
                path_name,
                record_start,
                seq_size,
                aln_is_reverse[seq_rank],
                path_length,
                aligned_seq
            });

            clear_string(path_name);
            clear_string(aligned_seq);
        }
    }

    // free memory
    free(is_rc);
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

    build_odgi_abPOA(ab, abpt, output_graph, names, aln_is_reverse, consensus_name, generate_consensus);

    abpoa_free(ab, abpt);
    abpoa_free_para(abpt);

    // normalize the representation, allowing for nodes > 1bp
    odgi::algorithms::unchop(output_graph);

    // order the graph
    output_graph.apply_ordering(odgi::algorithms::topological_order(&output_graph), true);

    // output_graph.to_gfa(std::cout);
    return output_graph;
}

odgi::graph_t smooth_spoa(const xg::XG &graph, const block_t &block,
                          const uint64_t &block_id,
                          std::int8_t poa_m, std::int8_t poa_n, std::int8_t poa_g,
                          std::int8_t poa_e, std::int8_t poa_q, std::int8_t poa_c,
                          bool local_alignment,
                          std::vector<maf_row_t> *maf, bool keep_sequence,
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

    std::string consensus;
    if (!consensus_name.empty()){
        consensus = poa_graph->generate_consensus();
        aln_is_reverse.push_back(false);  // the consensus is considered in forward
    }

    if (maf != nullptr) {
        std::vector<std::string> msa;
        poa_graph->generate_multiple_sequence_alignment(msa, !consensus_name.empty());

        uint64_t num_seqs = msa.size();
        for(uint64_t seq_rank = 0; seq_rank < num_seqs; seq_rank++){
            std::string path_name;
            uint64_t seq_size;
            uint64_t path_length;

            uint64_t record_start;
            if (consensus_name.empty() || seq_rank < num_seqs - 1){
                auto path_handle = graph.get_path_handle_of_step(block.path_ranges[seq_rank].begin);

                path_name = graph.get_path_name(path_handle);
                path_length = graph.get_path_length(path_handle);

                // If the strand field is "-" then this is the start relative to the reverse-complemented source sequence
                uint64_t path_range_begin = graph.get_position_of_step(block.path_ranges[seq_rank].begin);
                auto last_step = graph.get_previous_step(block.path_ranges[seq_rank].end);
                record_start = aln_is_reverse[seq_rank] ?
                               (path_length - graph.get_position_of_step(last_step) - graph.get_length(graph.get_handle_of_step(last_step))):
                               path_range_begin;

                seq_size = seqs[seq_rank].size(); // <==> block.path_ranges[seq_rank].length
            }else{
                // The last sequence is the gapped consensus

                path_name = consensus_name;
                path_length = consensus.size();
                record_start = 0;
                seq_size = path_length;
            }

            maf->push_back({
                path_name,
                record_start,
                seq_size,
                aln_is_reverse[seq_rank],
                path_length,
                keep_sequence ? msa[seq_rank] : ""
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
    build_odgi(poa_graph, output_graph, names, aln_is_reverse, consensus_name, !consensus_name.empty());

    // normalize the representation, allowing for nodes > 1bp
    odgi::algorithms::unchop(output_graph);

    // order the graph
    output_graph.apply_ordering(odgi::algorithms::topological_order(&output_graph), true);

    // output_graph.to_gfa(out);
    return output_graph;
}

void _put_block_in_group(
        maf_t &merged_maf_blocks, uint64_t block_id, uint64_t num_seq_in_block,
        bool add_consensus,
        std::vector<std::vector<maf_row_t>> &mafs
){
    merged_maf_blocks.field_blocks.push_back(block_id);

    for (uint64_t i = 0; i < num_seq_in_block; i ++){
        auto maf_row = mafs[block_id][i];

        if (!add_consensus || i < num_seq_in_block - 1) {
            merged_maf_blocks.rows.insert(std::pair<std::string, maf_partial_row_t>(
                    maf_row.path_name,
                    {
                            maf_row.record_start,
                            maf_row.seq_size,
                            maf_row.is_reversed,
                            maf_row.path_length,
                            maf_row.aligned_seq
                    }
            ));
        } else {
            merged_maf_blocks.consensus_rows.push_back(std::pair<std::string, maf_partial_row_t>(
                    maf_row.path_name,
                    {
                            maf_row.record_start,
                            maf_row.seq_size,
                            maf_row.is_reversed,
                            maf_row.path_length,
                            maf_row.aligned_seq
                    }
            ));
        }
    }
}

odgi::graph_t smooth_and_lace(const xg::XG &graph,
                              const std::vector<block_t> &blocks,
                              int poa_m, int poa_n,
                              int poa_g, int poa_e,
                              int poa_q, int poa_c,
                              bool local_alignment,
                              std::string &path_output_maf, std::string &maf_header, bool merge_blocks,
                              bool use_abpoa,
                              const std::string &consensus_base_name) {

    bool produce_maf = !path_output_maf.empty();
    bool add_consensus = !consensus_base_name.empty();

    //
    // record the start and end points of all the path ranges and the consensus
    //
    std::vector<odgi::graph_t> block_graphs;
    block_graphs.resize(blocks.size());

    // ToDo: better a queue (or nthreads queues)?
    // If merged consensus sequences have to be embedded, this structures are needed to keep the blocks' coordinates,
    // but the sequences will be considered (and kept in memory) only if a MAF has to be produced
    std::vector<std::vector<maf_row_t>> mafs(produce_maf || (add_consensus && merge_blocks) ? blocks.size() : 0);
    std::vector<std::atomic<bool>> mafs_ready(produce_maf || (add_consensus && merge_blocks) ? blocks.size() : 0);

    std::vector<std::pair<uint64_t, uint64_t>> merged_block_id_intervals;

    auto write_maf_lambda = [&]() {
        if (produce_maf || (add_consensus && merge_blocks)){
            uint64_t block_id = 0;

            maf_t merged_maf_blocks;

            uint64_t num_block = blocks.size();

            std::ofstream out_maf;

            if (produce_maf){
                out_maf.open(path_output_maf.c_str());
                out_maf << maf_header << std::endl;
            }

            while (block_id < num_block) {
                if (mafs_ready[block_id].load()){
                    uint64_t num_seq_in_block = mafs[block_id].size();
                    bool is_last_block = block_id == num_block - 1;

                    bool contains_loops = false;
                    std::unordered_set<path_handle_t> seen_paths;
                    for (auto &path_range : blocks[block_id].path_ranges){
                        path_handle_t path = graph.get_path_handle_of_step(path_range.begin);
                        if (seen_paths.count(path)) {
                            contains_loops = true;
                            break;
                        } else {
                            seen_paths.insert(path);
                        }
                    }
                    seen_paths.clear();

                    bool prep_new_merge_group = false;
                    bool merged = false;
                    if (merge_blocks){
                        if (!contains_loops) {
                            if (merged_maf_blocks.field_blocks.empty()){
                                merged = true;
                            } else {
                                // The merged group has to have the number of rows (then paths in this case) as the block to merge
                                if (merged_maf_blocks.rows.size() + (add_consensus ? 1 : 0) == num_seq_in_block) {
                                    merged = true;

                                    for (uint64_t i = 0; i < num_seq_in_block; i++) {
                                        // Do not check the consensus (always forward)
                                        if (!add_consensus || i < num_seq_in_block - 1) {
                                            auto maf_row = mafs[block_id][i];

                                            if (merged_maf_blocks.rows.count(maf_row.path_name) != 1){
                                                merged = false; // Current block not mergeable, so write the blocks which are waiting in memory
                                                prep_new_merge_group = !is_last_block;
                                                break;
                                            }

                                            maf_partial_row_t merged_maf_row = merged_maf_blocks.rows[maf_row.path_name];

                                            if (
                                                    (merged_maf_row.is_reversed != maf_row.is_reversed) ||
                                                    (merged_maf_row.is_reversed && ((merged_maf_row.record_start - maf_row.seq_size) != maf_row.record_start)) ||
                                                    (!merged_maf_row.is_reversed && ((merged_maf_row.record_start + merged_maf_row.seq_size) != maf_row.record_start))
                                                    ) {
                                                merged = false; // Current block not mergeable, so write the blocks which are waiting in memory
                                                prep_new_merge_group = !is_last_block;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (merged){
                        if (merged_maf_blocks.field_blocks.empty()){
                            _put_block_in_group(merged_maf_blocks, block_id, num_seq_in_block, add_consensus, mafs);
                        } else {
                            merged_maf_blocks.field_blocks.push_back(block_id);

                            for (uint64_t i = 0; i < num_seq_in_block; i ++){
                                auto maf_row = mafs[block_id][i];

                                if (!add_consensus || i < num_seq_in_block - 1) {
                                    if (maf_row.is_reversed){
                                        merged_maf_blocks.rows[maf_row.path_name].record_start -= maf_row.seq_size;

                                        merged_maf_blocks.rows[maf_row.path_name].aligned_seq = maf_row.aligned_seq + merged_maf_blocks.rows[maf_row.path_name].aligned_seq;
                                    }else{
                                        merged_maf_blocks.rows[maf_row.path_name].aligned_seq += maf_row.aligned_seq;
                                    }
                                    merged_maf_blocks.rows[maf_row.path_name].seq_size += maf_row.seq_size;
                                } else {
                                    // The merged consensus is created when the merged block is written into a file
                                    merged_maf_blocks.consensus_rows.push_back(std::pair<std::string, maf_partial_row_t>(
                                            maf_row.path_name,
                                            {
                                                    maf_row.record_start,
                                                    maf_row.seq_size,
                                                    maf_row.is_reversed,
                                                    maf_row.path_length,
                                                    maf_row.aligned_seq
                                            }
                                    ));
                                }
                            }
                        }

                        _clear_maf_block(block_id, mafs);
                    }

                    if (!merged || is_last_block) {
                        uint64_t merged_maf_blocks_size = merged_maf_blocks.field_blocks.size();

                        if (merged_maf_blocks_size > 0){
                            /*std::stringstream joined_block_ids;
                            for (size_t i = 0; i < merged_maf_blocks_size; ++i) {
                                if (i != 0) { joined_block_ids << ","; }
                                joined_block_ids << merged_maf_blocks.field_blocks[i];
                            }*/

                            std::string block_id_range = std::to_string(merged_maf_blocks.field_blocks[0]);
                            if (merged_maf_blocks_size > 1) {
                                block_id_range += "-" + std::to_string(merged_maf_blocks.field_blocks[merged_maf_blocks_size - 1]);

                                merged_block_id_intervals.emplace_back(
                                        merged_maf_blocks.field_blocks[0], merged_maf_blocks.field_blocks[merged_maf_blocks_size - 1]
                                );
                            }

                            if (produce_maf){
                                out_maf << "a blocks=" << block_id_range << " loops=false";
                                if (merged_maf_blocks_size > 1) {
                                    out_maf << " merged=true";
                                }
                                out_maf << std::endl;

                                std::vector<maf_row_t> rows;
                                for (auto& maf_row : merged_maf_blocks.rows) {
                                    rows.push_back(
                                            {
                                                    maf_row.first,
                                                    maf_row.second.record_start,
                                                    maf_row.second.seq_size,
                                                    maf_row.second.is_reversed,
                                                    maf_row.second.path_length,
                                                    maf_row.second.aligned_seq
                                            }
                                    );
                                }
                                if (add_consensus){
                                    uint64_t merged_consensus_seq_size = 0;
                                    uint64_t merged_consensus_path_length = 0;
                                    std::string merged_consensus_aligned_seq;

                                    uint64_t length_alignment = merged_maf_blocks.rows.begin()->second.aligned_seq.length();
                                    uint64_t start_cons_pos_in_alignment = 0;
                                    for (auto & maf_cons_row : merged_maf_blocks.consensus_rows){
                                        std::string gapped_cons;

                                        uint64_t pos;
                                        for (pos = 0; pos < start_cons_pos_in_alignment; pos++){ gapped_cons += "-"; }

                                        gapped_cons += maf_cons_row.second.aligned_seq;

                                        pos += maf_cons_row.second.aligned_seq.length();
                                        for (; pos < length_alignment; pos++){ gapped_cons += "-"; }

                                        rows.push_back(
                                                {
                                                        maf_cons_row.first,
                                                        maf_cons_row.second.record_start,
                                                        maf_cons_row.second.seq_size,
                                                        maf_cons_row.second.is_reversed,
                                                        maf_cons_row.second.path_length,
                                                        gapped_cons
                                                }
                                        );

                                        clear_string(gapped_cons);

                                        start_cons_pos_in_alignment += maf_cons_row.second.aligned_seq.length();


                                        if (merged_maf_blocks_size > 1) {
                                            merged_consensus_seq_size += maf_cons_row.second.seq_size;
                                            merged_consensus_path_length += maf_cons_row.second.path_length;
                                            merged_consensus_aligned_seq += maf_cons_row.second.aligned_seq;
                                        }
                                    }

                                    if (merged_maf_blocks_size > 1) {
                                        // Write the merged consensus
                                        rows.push_back(
                                                {
                                                        consensus_base_name + block_id_range + " ",
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
                                write_maf_rows(out_maf, rows);

                                for(auto& row : rows){
                                    clear_string(row.path_name);
                                    clear_string(row.aligned_seq);
                                }
                                clear_vector(rows);
                            }

                            clear_string(block_id_range);

                            clear_vector(merged_maf_blocks.field_blocks);
                            for (auto& maf_row : merged_maf_blocks.rows){
                                clear_string(maf_row.first);
                                clear_string(maf_row.second.aligned_seq);
                            }
                            merged_maf_blocks.rows.clear();

                            for (auto& maf_cons_row : merged_maf_blocks.consensus_rows){
                                clear_string(maf_cons_row.first);
                                clear_string(maf_cons_row.second.aligned_seq);
                            }
                            clear_vector(merged_maf_blocks.consensus_rows);
                        }

                        if (!merged && !prep_new_merge_group) {
                            if (produce_maf){
                                out_maf << "a blocks=" + std::to_string(block_id) << " loops=" << (contains_loops ? "true" : "false") << std::endl;
                                write_maf_rows(out_maf, mafs[block_id]);
                            }

                            _clear_maf_block(block_id, mafs);
                        }
                    }

                    if (prep_new_merge_group){
                        // This is a mergeable (and not the last) block: it is the first one, or the last merge failed
                        // (and the current un-merged block becomes the first one of the next group)

                        _put_block_in_group(merged_maf_blocks, block_id, num_seq_in_block, add_consensus, mafs);

                        _clear_maf_block(block_id, mafs);
                    }

                    block_id++;
                }

                std::this_thread::sleep_for(std::chrono::nanoseconds(1));
            }

            if (produce_maf){
                out_maf.close();
            }

            clear_vector(mafs);
            clear_vector(mafs_ready);
        }
    };
    std::thread write_maf_thread(write_maf_lambda);

    std::vector<path_position_range_t> path_mapping;
    std::vector<path_position_range_t> consensus_mapping;
    std::mutex path_mapping_mutex, consensus_mapping_mutex, logging_mutex;
    uint64_t thread_count = odgi::get_thread_count();

    std::stringstream poa_banner;
    poa_banner << "[smoothxg::smooth_and_lace] applying "
               << (local_alignment ? "local" : "global") << " "
               << (use_abpoa ? "abPOA" : "SPOA")
               << " to " << blocks.size() << " blocks:";
    progress_meter::ProgressMeter poa_progress(blocks.size(), poa_banner.str());

    paryfor::parallel_for<uint64_t>(
        0, blocks.size(), thread_count, [&](uint64_t block_id, int tid) {
            auto &block = blocks[block_id];

            std::string consensus_name;
            if (add_consensus){
                consensus_name = consensus_base_name + std::to_string(block_id);
            }

            // std::cerr << "on block " << block_id+1 << " of " << blocks.size() << std::endl;
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
                                           (produce_maf || (add_consensus && merge_blocks)) ? &mafs[block_id] : nullptr,
                                           produce_maf,
                                           true, // banded alignment
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
                                          (produce_maf || (add_consensus && merge_blocks)) ? &mafs[block_id] : nullptr,
                                          produce_maf,
                                          consensus_name);
            }

            if (produce_maf || (add_consensus && merge_blocks)){
                mafs_ready[block_id].store(true);
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
                    auto consensus_handle = block_graph.get_path_handle(consensus_name);

                    uint64_t path_end = 0;
                    step_handle_t empty_step;
                    as_integers(empty_step)[0] = 0;
                    as_integers(empty_step)[1] = 0;
                    block_graph.for_each_step_in_path(consensus_handle, [&](const step_handle_t &step) {
                        path_end += block_graph.get_length(block_graph.get_handle_of_step(step));
                    });

                    {
                        std::lock_guard<std::mutex> guard(consensus_mapping_mutex);
                        consensus_mapping.push_back(
                                {
                                    as_path_handle(0),  // consensus = 0 path handle
                                    0,                        // start position
                                    path_end,                          // end position
                                    empty_step, empty_step, consensus_handle,
                                    block_id
                                }
                        );
                    }
                }
            }
            poa_progress.increment(1);
        });

    poa_progress.finish();

    write_maf_thread.join();

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

    std::stringstream add_graph_banner;
    add_graph_banner << "[smoothxg::smooth_and_lace] adding "
                     << block_graphs.size() << " graphs:";
    progress_meter::ProgressMeter add_graph_progress(block_graphs.size(), add_graph_banner.str());

    for (auto &block : block_graphs) {
        uint64_t id_trans = smoothed.get_node_count();
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
        add_graph_progress.increment(1);
    }
    add_graph_progress.finish();
    
    // then for each path, ensure that it's embedded in the graph by walking
    // through its block segments in order and linking them up in the output
    // graph
    std::stringstream lace_banner;
    lace_banner << "[smoothxg::smooth_and_lace] embedding "
                << path_mapping.size() << " path fragments:";
    progress_meter::ProgressMeter lace_progress(path_mapping.size(), lace_banner.str());
    for (uint64_t i = 0; i < path_mapping.size(); ++i) {
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
            lace_progress.increment(1);
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
    lace_progress.finish();
    // now verify that smoothed has paths that are equal to the base graph
    // and that all the paths are fully embedded in the graph
    std::vector<path_handle_t> paths;
    smoothed.for_each_path_handle(
        [&](const path_handle_t &path) {
            paths.push_back(path);
        });
    std::stringstream validate_banner;
    validate_banner << "[smoothxg::smooth_and_lace] validating "
                    << paths.size() << " path sequences:";
    progress_meter::ProgressMeter validate_progress(paths.size(), validate_banner.str());
    paryfor::parallel_for<uint64_t>(
        0, paths.size(), thread_count,
        [&](uint64_t path_id, int tid) {
            auto path = paths[path_id];
            std::string orig_seq, smoothed_seq;
            graph.for_each_step_in_path(
                graph.get_path_handle(smoothed.get_path_name(path)),
                [&](const step_handle_t &step) {
                    orig_seq.append(
                        graph.get_sequence(graph.get_handle_of_step(step)));
                });
            smoothed.for_each_step_in_path(
                path,
                [&](const step_handle_t &step) {
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
            validate_progress.increment(1);
        });
    validate_progress.finish();

    if (!consensus_mapping.empty()) {
        std::cerr << "[smoothxg::smooth_and_lace] sorting consensus"
                  << std::endl;

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

        std::cerr << "[smoothxg::smooth_and_lace] embedding consensus"
                  << std::endl;

        for (auto &pos_range : consensus_mapping) {
            auto &block = block_graphs[pos_range.target_graph_id];
            path_handle_t smoothed_path = smoothed.create_path_handle(
                block.get_path_name(pos_range.target_path)
            );

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

        // Sort respect to the target_graph_id (== block_id), because in the merged_block_id_intervals
        // there are the block_id_intervals expressed as first_block_id and last_block_id
        ips4o::parallel::sort(consensus_mapping.begin(), consensus_mapping.end(),
                              [](const path_position_range_t &a, const path_position_range_t &b) {
                                  return (a.target_graph_id < b.target_graph_id);
                              });
        for (auto& merged_block_id_interval : merged_block_id_intervals){
            path_handle_t consensus_path = smoothed.create_path_handle(
                    consensus_base_name +
                    std::to_string(merged_block_id_interval.first) + "-" + std::to_string(merged_block_id_interval.second)
            );

            for (uint64_t block_id = merged_block_id_interval.first; block_id <= merged_block_id_interval.second; block_id++) {
                auto &block = block_graphs[block_id];
                auto &id_trans = id_mapping[block_id];
                block.for_each_step_in_path(consensus_mapping[block_id].target_path, [&](const step_handle_t &step) {
                    handle_t h = block.get_handle_of_step(step);
                    handle_t t = smoothed.get_handle(block.get_id(h) + id_trans, block.get_is_reverse(h));
                    smoothed.append_step(consensus_path, t);
                });
            }
        }

        // todo: validate the consensus paths as well
    }

    std::stringstream embed_banner;
    embed_banner << "[smoothxg::smooth_and_lace] walking edges in "
                 << paths.size() << " paths:";
    progress_meter::ProgressMeter embed_progress(paths.size(), embed_banner.str());
    // embed all paths in the graph
    smoothed.for_each_path_handle(
        [&](const path_handle_t& path) {
            handle_t last;
            step_handle_t begin_step = smoothed.path_begin(path);
            smoothed.for_each_step_in_path(
                path,
                [&](const step_handle_t &step) {
                    handle_t h = smoothed.get_handle_of_step(step);
                    if (step != begin_step) {
                        smoothed.create_edge(last, h);
                    }
                    last = h;
                });
            embed_progress.increment(1);
        });
    embed_progress.finish();

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
