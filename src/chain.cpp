#include "chain.hpp"

namespace smoothxg {

std::vector<anchor_t> anchors_for_path(
    const xg::XG& graph,
    const path_handle_t& path) {
    std::vector<anchor_t> anchors;
    graph.for_each_step_in_path(
        path,
        [&](const step_handle_t& query_step) {
            seq_pos_t query_position = seq_pos::encode(graph.get_position_of_step(query_step), false);
            const handle_t& handle = graph.get_handle_of_step(query_step);
            size_t handle_length = graph.get_length(handle);
            //size_t query_is_rev = graph.get_is_reverse(handle);
            graph.for_each_step_position_on_handle(
                handle,
                [&](const step_handle_t& target_step,
                    const bool& target_rev_vs_query,
                    const size_t& target_pos) {
                    path_handle_t target_path = graph.get_path_handle_of_step(target_step);
                    seq_pos_t target_position = seq_pos::encode(target_pos, target_rev_vs_query);
                    if (query_step != target_step) {
                        anchors.emplace_back(
                            path,
                            query_position,
                            query_position + handle_length,
                            target_path,
                            target_position,
                            target_position + handle_length);
                    }
                    return true;
                });
        });
    return anchors;
}

/*
std::vector<anchor_t>
anchors_for_query(
    const gyeet_index_t& index,
    const char* seq,
    const size_t& len) {
    std::vector<anchor_t> anchors;
    // for each kmer
    const uint64_t& kmer_length = index.kmer_length;
    for (uint64_t i = 0; i <= len-kmer_length; ++i) {
        //std::cerr << std::string(seq+i, kmer_length) << " " << i << " " << std::endl;
        index.for_values_of(
            seq + i, kmer_length,
            [&](const kmer_start_end_t& v) {
                std::cerr << std::string(seq+i, kmer_length) << " " << i << " "
                          << seq_pos::offset(v.begin) << " " << seq_pos::offset(v.end)
                          << std::endl;
                anchors.emplace_back(seq_pos::encode(i, false),
                                     seq_pos::encode(i+kmer_length, false),
                                     v.begin, v.end);
            });
    }
    return anchors;
}
*/

void chain_t::compute_boundaries(const double& mismatch_rate) {
    // find the longest contiguous range covered by the chain
    // where the size is no more than some function of our seed_length * number of anchors * 1+mismatch_rate
    /*
    for (auto& anchor : anchors) {
        std::cerr << "\"" << &anchor << "\" "
                  << "[shape=box label=\"" << "score=" << anchor->max_chain_score << ","
                  << "(" << seq_pos::to_string(anchor->query_begin) << ".."
                  << seq_pos::to_string(anchor->query_end) << "),"
                  << "(" << seq_pos::to_string(anchor->target_begin) << ".."
                  << seq_pos::to_string(anchor->target_end) << ")\"];" << std::endl;
    }
    */
    // assign target path
    target_path = anchors.front()->target_path;
    // find the inner boundaries
    const seq_pos_t& first_target_end = anchors.front()->target_end;
    const seq_pos_t& last_target_begin = anchors.back()->target_begin;
    // can we linearly extend to the start, end?

    const seq_pos_t& first_target_begin = anchors.front()->target_begin;
    const seq_pos_t& last_target_end = anchors.back()->target_end;
    //if (first_target_begin < last_target_end) {
    if (seq_pos::is_rev(first_target_begin) == seq_pos::is_rev(last_target_end) &&
        first_target_begin < last_target_end && score * (1+mismatch_rate) > last_target_end - first_target_begin) {
        target_begin = first_target_begin;
        target_end = last_target_end;
    } else if (seq_pos::is_rev(first_target_end) == seq_pos::is_rev(last_target_begin)
               && first_target_end < last_target_begin) {
        target_begin = first_target_end;
        target_end = last_target_begin;
    } else {
        // kill chains that jump nonlinearly internally
        //target_begin = first_target_end;
        //target_end = last_target_begin;
        score = -std::numeric_limits<double>::max();
    }
    //std::cerr << "target begin " << seq_pos::to_string(target_begin) << std::endl;
    //std::cerr << "target end " << seq_pos::to_string(target_end) << std::endl;
}

std::vector<chain_t>
chains(const std::string& query_name,
       std::vector<anchor_t>& anchors,
       const uint64_t& max_gap,
       const double& mismatch_rate,
       const uint64_t& chain_min_n_anchors,
       const uint64_t bandwidth,
       const double secondary_chain_threshold,
       const double max_mapq) {
    // sort the anchors by their ending position in the target
    std::sort(anchors.begin(),
              anchors.end(),
              [](const anchor_t& a,
                 const anchor_t& b) {
                  return as_integer(a.target_path) < as_integer(b.target_path)
                        || (a.target_path == b.target_path && a.target_end < b.target_end);
              });
    // calculate max chaining scores under our bandwidth
    // we walk from the last to first anchor
    for (int64_t i = 0; i < anchors.size(); ++i) {
        anchor_t& anchor_i = anchors[i];
        //std::cerr << "anchor_i " << anchor_i.query_begin << ".." << anchor_i.query_end << std::endl;
        anchor_i.max_chain_score = 1; // XXXXXXX TODO // seed_length;
        int64_t min_j = bandwidth > i ? 0 : i - bandwidth;
        for (int64_t j = i-1; j >= min_j; --j) {
            anchor_t& anchor_j = anchors[j];
            //std::cerr << "anchor_j " << anchor_j.query_begin << ".." << anchor_j.query_end << std::endl;
            //if (anchor_i.target_end - anchor_j.target_end > max_gap) break;
            double proposed_score = score_anchors(anchor_j,
                                                  anchor_i,
                                                  max_gap);
            //std::cerr << "proposed " << proposed_score << " vs " << anchor_i.max_chain_score << std::endl;
            //std::cerr << "diff " << proposed_score - anchor_i.max_chain_score << std::endl;
            if (proposed_score > anchor_i.max_chain_score) {
                //std::cerr << "taken!" << std::endl;
                anchor_i.max_chain_score = proposed_score;
                anchor_i.best_predecessor = &anchor_j;
            }
        }
    }

    std::ofstream out(query_name + ".chains.dot");
    out << "digraph G {" << std::endl;
    //out << "rankdir=LR;" << std::endl;
    for (auto& anchor : anchors) {
        out << "\"" << &anchor << "\" "
            << "[shape=box label=\"" << "score=" << anchor.max_chain_score << ","
            << "(" << seq_pos::to_string(anchor.query_begin) << ".."
            << seq_pos::to_string(anchor.query_end) << "),"
            << "(" << seq_pos::to_string(anchor.target_begin) << ".."
            << seq_pos::to_string(anchor.target_end) << ")\"];" << std::endl
            << "\"" << anchor.best_predecessor << "\" -> \"" << &anchor << "\";" << std::endl;
    }
    out << "}" << std::endl;
    out.close();

    // collect chains
    std::vector<chain_t> chains;
    int64_t i = anchors.size()-1;
    while (i >= 0) {
        anchor_t* a = &anchors[i];
        //std::cerr << "best predecessor " << a->best_predecessor << " "<< a->max_chain_score <<std::endl;
        if (a->best_predecessor != nullptr
            && a->max_chain_score > 1) { /// XXXXXXX FIXME // seed_length) { //!curr_chain) {
            chains.emplace_back();
            auto& curr_chain = chains.back();
            curr_chain.anchors.push_back(a);
            curr_chain.score = a->max_chain_score;
            do {
                anchor_t* b = a->best_predecessor;
                curr_chain.anchors.push_back(b);
                a->best_predecessor = nullptr; // mark done
                a = b; // swap
            } while (a->best_predecessor != nullptr);
            // investigate if a deque would be faster
            if (curr_chain.anchors.size() < chain_min_n_anchors) {
                chains.pop_back();
            } else {
                std::reverse(curr_chain.anchors.begin(),
                             curr_chain.anchors.end());
            }
        }
        --i;
    }
    //std::cerr << "chain count " <<chains.size() <<std::endl;
    // sort the chains by score, descending
    std::sort(chains.begin(), chains.end(),
              [](const chain_t& a,
                 const chain_t& b) {
                  return a.score > b.score;
              });
    /*
    uint64_t j = 0;
    for (auto& chain : chains) {
        std::cerr << "chain " << ++j << " " << &chain << " " << chain.anchors.size() << " " << chain.score << " " << chain.mapping_quality << " " << chain.is_secondary << " " << chain.processed() << std::endl;
    }
    */
    // find the primary chains by examining their overlaps in the query
    IITree<seq_pos_t, chain_t*> tree;
    for (auto& chain : chains) {
        const seq_pos_t& query_begin = chain.anchors.front()->query_begin;
        const seq_pos_t& query_end = chain.anchors.back()->query_end;
        tree.add(query_begin, query_end, &chain);
    }
    tree.index();
    for (auto& chain : chains) {
        if (!chain.processed()) {
            const seq_pos_t& chain_begin = chain.anchors.front()->query_begin;
            const seq_pos_t& chain_end = chain.anchors.back()->query_end;
            uint64_t chain_length = chain_end - chain_begin;
            chain_t* best_secondary = nullptr;
            std::vector<size_t> ovlp;
            tree.overlap(chain_begin, chain_end, ovlp);
            for (auto& idx : ovlp) {
                chain_t* other_chain = tree.data(idx);
                if (other_chain->target_path == chain.target_path
                    && other_chain != &chain && other_chain->score <= chain.score) {
                    const seq_pos_t& other_begin = tree.start(idx);
                    const seq_pos_t& other_end = tree.end(idx);
                    uint64_t other_length = other_end - other_begin;
                    seq_pos_t ovlp_begin = std::max(chain_begin, other_begin);
                    seq_pos_t ovlp_end = std::min(chain_end, other_end);
                    uint64_t ovlp_length = ovlp_end - ovlp_begin;
                    if ((double)ovlp_length > (double)other_length * secondary_chain_threshold) {
                        // this chain is secondary
                        other_chain->mapping_quality = 0;
                        other_chain->is_secondary = true;
                    }
                    if (best_secondary == nullptr
                        || best_secondary->score < other_chain->score) {
                        best_secondary = other_chain;
                    }
                }
            }
            if (best_secondary == nullptr) {
                chain.mapping_quality = max_mapq;
            } else {
                chain.mapping_quality =
                    40 * (1 - best_secondary->score / chain.score)
                    * std::min(1.0, (double)chain.anchors.size()/10.0)
                    * log(chain.score);
            }
        }
    }
    for (auto& chain : chains) {
        chain.compute_boundaries(mismatch_rate);
    }
    return chains;
}

double score_anchors(const anchor_t& a,
                     const anchor_t& b,
                     const uint64_t& max_gap) {
    if (a.query_path != b.query_path ||
        a.target_path != b.target_path ||
        a.query_end >= b.query_end ||
        !( seq_pos::is_rev(a.target_end)
           == seq_pos::is_rev(b.target_end)
           == seq_pos::is_rev(a.target_begin)
           == seq_pos::is_rev(b.target_begin))) {
        return -std::numeric_limits<double>::max();
    } else {
        uint64_t query_length = std::min(b.query_begin - a.query_begin,
                                         b.query_end - a.query_end);
        uint64_t query_overlap = (a.query_end > b.query_begin ? a.query_end - b.query_begin : 0);
        uint64_t target_length = std::min(b.target_begin - a.target_begin,
                                          b.target_end - a.target_end);
        //std::cerr << "query_length " << query_length << " target_length " << target_length << std::endl;
        uint64_t gap_length = std::abs((int64_t)query_length - (int64_t)target_length);
        if (gap_length > max_gap) { //std::max(query_length, target_length) > max_gap) {
            return -std::numeric_limits<double>::max();
            //return std::numeric_limits<double>::min();
        } else {
            //std::cerr << "query_length " << query_length << " target_length " << target_length << std::endl;
            double gap_cost = gap_length == 0 ? 0
                : 0.01 * gap_length + 0.5 * log2(gap_length);
                //: 0.01 * seed_length * gap_length + 0.5 * log2(gap_length);
            uint64_t match_length = std::min(query_length, target_length);
            //std::cerr << "chain score is " << a.max_chain_score + match_length - gap_cost << std::endl;
            // round to 3 decimal digits, avoid problems with floating point instability causing chain truncation
            return std::round((a.max_chain_score + match_length - gap_cost) * 1000.0) / 1000.0 + query_overlap;
        }
    }
}

uint64_t chain_query_length(const chain_t& chain) {
    return chain.anchors.back()->query_end - chain.anchors.front()->query_begin;
}

std::vector<superchain_t>
superchains(const std::string& query_name,
            std::vector<chain_t>& chains,
            const double& mismatch_rate,
            const double& chain_overlap_max,
            const uint64_t bandwidth) {
    // sort the chains by end coordinate in the query
    // score each with the previous N
    std::vector<chain_node_t> chain_nodes; //chain_nodes.reserve(chains.size());
    // collect primary chains
    for (auto& chain : chains) {
        if (!chain.is_secondary) {
            chain_nodes.emplace_back(&chain);
        }
    }
    // sort by end in query
    std::sort(chain_nodes.begin(), chain_nodes.end(),
              [](const chain_node_t& a,
                 const chain_node_t& b) {
                  return a.chain->anchors.back()->query_end
                      < b.chain->anchors.back()->query_end;
              });
    // dynamic programming
    for (int64_t i = 0; i < chain_nodes.size(); ++i) {
        chain_node_t& chain_node_i = chain_nodes[i];
        //std::cerr << "anchor_i " << anchor_i.query_begin << ".." << anchor_i.query_end << std::endl;
        chain_node_i.max_superchain_score = chain_node_i.chain->score; // chain_query_length(*chain_node_i.chain)
        int64_t min_j = bandwidth > i ? 0 : i - bandwidth;
        for (int64_t j = i-1; j >= min_j; --j) {
            chain_node_t& chain_node_j = chain_nodes[j];
            //anchor_t& anchor_j = anchors[j];
            //std::cerr << "anchor_j " << anchor_j.query_begin << ".." << anchor_j.query_end << std::endl;
            double proposed_score = score_chain_nodes(chain_node_j, chain_node_i, chain_overlap_max);
            //std::cerr << "proposed " << proposed_score << " vs " << chain_node_i.max_superchain_score << std::endl;
            //std::cerr << "diff " << proposed_score - chain_node_i.max_superchain_score << std::endl;
            if (proposed_score > chain_node_i.max_superchain_score) {
                //std::cerr << "taken!" << std::endl;
                chain_node_i.max_superchain_score = proposed_score;
                chain_node_i.best_predecessor = &chain_node_j;
            }
        }
    }

    std::ofstream out(query_name + ".superchains.dot");
    out << "digraph G {" << std::endl;
    //out << "rankdir=LR;" << std::endl;
    for (auto& chain_node : chain_nodes) {
        auto& first_anchor = *chain_node.chain->anchors.front();
        auto& last_anchor = *chain_node.chain->anchors.back();
        out << "\"" << &chain_node << "\" "
            << "[shape=box label=\""
            << "anchors=" << chain_node.chain->anchors.size() << ","
            << "score=" << chain_node.max_superchain_score << ","
            << "(" << seq_pos::to_string(first_anchor.query_begin) << ".."
            << seq_pos::to_string(last_anchor.query_end) << "),"
            << "(" << seq_pos::to_string(first_anchor.target_begin) << ".."
            << seq_pos::to_string(last_anchor.target_end) << ")\"];" << std::endl
            << "\"" << chain_node.best_predecessor << "\" -> \"" << &chain_node << "\";" << std::endl;
    }
    out << "}" << std::endl;
    out.close();

    // collect superchains
    std::vector<superchain_t> superchains;
    int64_t i = chain_nodes.size()-1;
    while (i >= 0) {
        chain_node_t* a = &chain_nodes[i];
        //std::cerr << "best predecessor " << a->best_predecessor << " "<< a->max_chain_score <<std::endl;
        if (a->best_predecessor != nullptr
            && a->max_superchain_score) { // HMMMM
            //std::cerr << "adding superchain" << std::endl;
            superchains.emplace_back();
            auto& curr_superchain = superchains.back();
            curr_superchain.is_secondary = curr_superchain.is_secondary && a->chain->is_secondary;
            curr_superchain.chains.push_back(a->chain);
            a->used = true;
            curr_superchain.score = a->max_superchain_score;
            do {
                chain_node_t* b = a->best_predecessor;
                curr_superchain.chains.push_back(b->chain);
                curr_superchain.is_secondary = curr_superchain.is_secondary && b->chain->is_secondary;
                b->used = true;
                a->best_predecessor = nullptr; // mark done
                a = b; // swap
            } while (a->best_predecessor != nullptr);
            // investigate if a deque would be faster
            std::reverse(curr_superchain.chains.begin(),
                         curr_superchain.chains.end());
        }
        --i;
    }
    // add back in all chains that aren't part of a superchain
    for (auto& chain_node : chain_nodes) {
        if (!chain_node.used) {
            superchains.emplace_back();
            auto& superchain = superchains.back();
            superchain.is_secondary = chain_node.chain->is_secondary;
            superchain.chains.push_back(chain_node.chain);
            superchain.score = chain_node.chain->score;
        }
    }
    // re-score the chains
    for (auto& superchain : superchains) {
        superchain.score = 0;
        for (uint64_t i = 0; i < superchain.chains.size()-1; ++i) {
            superchain.score = std::max(superchain.score,
                                        score_chain_pair(*superchain.chains[i],
                                                         superchain.score,
                                                         *superchain.chains[i+1],
                                                         superchain.chains[i+1]->score,
                                                         chain_overlap_max));
        }
    }
    std::cerr << "superchains for " << query_name << std::endl;
    std::cerr << "chain count " << chains.size() <<std::endl;
    std::cerr << "superchain count " << superchains.size() <<std::endl;
    // sort the superchains by score, descending
    std::sort(superchains.begin(), superchains.end(),
              [](const superchain_t& a,
                 const superchain_t& b) {
                  return a.score > b.score;
              });
    return superchains;
}

double score_chain_pair(const chain_t& a,
                        const double& score_a,
                        const chain_t& b,
                        const double& score_b,
                        const double& overlap_max) {
    // ignore if the a chain is contained in b by more than our overlap_max
    if (a.target_path != b.target_path
        || (int64_t)seq_pos::offset(a.query_end()) - (int64_t)seq_pos::offset(b.query_begin())
        > std::round((b.query_end() - b.query_begin()) * overlap_max)) {
        return -std::numeric_limits<double>::max();
    } else {
        // otherwise, add scores
        return score_a + score_b;
    }
}
    
double score_chain_nodes(const chain_node_t& a,
                         const chain_node_t& b,
                         const double& overlap_max) {
    return score_chain_pair(*a.chain, a.max_superchain_score, *b.chain, b.chain->score, overlap_max);
}

/*
GAF format
https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf

Col     Type    Description
1       string  Query sequence name
2       int     Query sequence length
3       int     Query start (0-based; closed)
4       int     Query end (0-based; open)
5       char    Strand relative to the path: "+" or "-"
6       string  Path matching /([><][^\s><]+(:\d+-\d+)?)+|([^\s><]+)/
7       int     Path length
8       int     Start position on the path (0-based)
9       int     End position on the path (0-based)
10      int     Number of residue matches
11      int     Alignment block length
12      int     Mapping quality (0-255; 255 for missing)
*/

void for_handle_at_anchor_begin_in_chain(
    const chain_t& chain,
    const xg::XG& index,
    const std::function<void(const handle_t&)>& func) {
    handle_t last = max_handle();
    for (auto& anchor : chain.anchors) {
        handle_t curr = index.get_handle_of_step(
            index.get_step_at_position(anchor->target_path,
                                       seq_pos::offset(anchor->target_begin)));
        if (curr != last) {
            func(curr);
            last = curr;
        }
    }
}

void write_chain_gaf(
    std::ostream& out,
    const chain_t& chain,
    const xg::XG& index,
    const std::string& query_name,
    const uint64_t& query_length) {
    out << query_name << "\t"
        << query_length << "\t"
        << chain.anchors.front()->query_begin << "\t"
        << chain.anchors.back()->query_end << "\t"
        << "+" << "\t";  // we're always forward strand relative to our path
    uint64_t path_length = 0;
    for_handle_at_anchor_begin_in_chain(
        chain, index,
        [&out,&index,&path_length](const handle_t& h) {
            path_length += index.get_length(h);
            out << (handle_is_rev(h) ? "<" : ">") << to_id(h);
        });
    out << "\t"
        << path_length << "\t"
        << 0 << "\t" // fixme
        << path_length << "\t" // fixme
        << path_length << "\t" // fixme
        << path_length << "\t" // fixme
        << std::min((int)std::round(chain.mapping_quality), 254) << "\t"
        << "ta:Z:chain" << std::endl;
}

void write_superchain_gaf(
    std::ostream& out,
    const superchain_t& superchain,
    const xg::XG& index,
    const std::string& query_name,
    const uint64_t& query_length) {
    out << query_name << "\t"
        << query_length << "\t"
        << superchain.chains.front()->anchors.front()->query_begin << "\t"
        << superchain.chains.back()->anchors.back()->query_end << "\t"
        << "+" << "\t";  // we're always forward strand relative to our path
    uint64_t path_length = 0;
    for (auto& chain : superchain.chains) {
        for_handle_at_anchor_begin_in_chain(
            *chain, index,
            [&out,&index,&path_length](const handle_t& h) {
                path_length += index.get_length(h);
                out << (handle_is_rev(h) ? "<" : ">") << to_id(h);
            });
    }
    out << "\t"
        << path_length << "\t"
        << 0 << "\t" // fixme
        << path_length << "\t" // fixme
        << path_length << "\t" // fixme
        << path_length << "\t" // fixme
        << 0 << "\t"
        << "ta:Z:superchain" << std::endl;
        //<< std::min((int)std::round(superchain.mapping_quality), 254) << std::endl;
}

std::vector<std::pair<path_handle_t, std::vector<superchain_t>>>
collinear_blocks(
    const xg::XG& graph,
    const uint64_t& max_gap,
    const double& mismatch_rate,
    const uint64_t& chain_min_n_anchors,
    const double& chain_overlap_max) {

    std::vector<std::pair<path_handle_t, std::vector<superchain_t>>> blocks;
    graph.for_each_path_handle(
        [&](const path_handle_t& path) {
            auto anchors = anchors_for_path(graph, path);
            std::string path_name = graph.get_path_name(path);
            uint64_t path_length = graph.get_path_length(path);
            auto query_chains = chains(path_name,
                                       anchors,
                                       max_gap,
                                       mismatch_rate,
                                       chain_min_n_anchors);
            for (auto& chain : query_chains) {
                write_chain_gaf(std::cerr, chain, graph, path_name, path_length);
            }
            auto query_superchains = superchains(path_name,
                                                 query_chains,
                                                 mismatch_rate,
                                                 chain_overlap_max);
            for (auto& superchain : query_superchains) {
                std::cerr << path_name << " -- " << graph.get_path_name(superchain.chains.front()->target_path) << " length " << superchain.chains.size() << " score " << superchain.score << std::endl;
               }
            blocks.emplace_back(std::make_pair(path, query_superchains));
        });
    return blocks;
}

}
