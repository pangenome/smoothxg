#include <deps/sdsl-lite/include/sdsl/int_vector.hpp>
#include "blocks.hpp"
#include "progress.hpp"

namespace smoothxg {

void smoothable_blocks(
    const odgi::graph_t& graph,
    blockset_t& blockset,
    const uint64_t& max_block_weight,
    const uint64_t& max_block_path_length,
    const uint64_t& max_path_jump,
    const uint64_t& max_edge_jump,
    const bool& order_paths_from_longest,
    const int num_threads,
	const odgi::algorithms::step_index_t &step_index,
	const std::vector<uint64_t> &node_offsets
    ) {
    // iterate over the handles in their vectorized order, collecting blocks that we can potentially smooth
    block_t block;
    std::vector<handle_t> block_handles;
    std::vector<sdsl::bit_vector> seen_steps(graph.get_path_count());

    // cast to vectorizable graph for determining the sort position of nodes
    std::cerr << "[smoothxg::smoothable_blocks] computing blocks" << std::endl;

    uint64_t rank = 0;
    graph.for_each_path_handle([&](const path_handle_t& path) {
        sdsl::util::assign(seen_steps[rank++], sdsl::bit_vector(graph.get_step_count(path), 0));
    });

    rank = 0;

    auto seen_step =
        [&](const step_handle_t& step) {
            return seen_steps[as_integer(graph.get_path_handle_of_step(step)) - 1][step_index.get_position(step, graph)];
        };
    auto mark_step =
        [&](const step_handle_t& step) {
            seen_steps[as_integer(graph.get_path_handle_of_step(step)) - 1][step_index.get_position(step, graph)] = 1;
        };
    auto toposplit_block =
        [&](const block_t& block) {
            ska::flat_hash_map<uint64_t, uint64_t> id_to_entry;
            std::vector<uint64_t> nodes;
            // collect handles
            uint64_t node_count = 0;
            for (auto& path_range : block.path_ranges) {
				step_handle_t step = path_range.begin;
                while (true) {
                    auto h = graph.get_handle_of_step(step);
                    auto id = graph.get_id(h);
                    auto f = id_to_entry.find(id);
                    if (f == id_to_entry.end()) {
                        id_to_entry[id] = node_count++;
                        nodes.push_back(id);
                    }
					if (step == path_range.end) break;
					step = graph.get_next_step(step);
                }
            }
            std::vector<std::atomic<odgi::DisjointSets::Aint>> data(nodes.size()+1); // maps into this set of disjoint sets
            auto dset = odgi::DisjointSets(data.data(), data.size());
            for (auto& path_range : block.path_ranges) {
                step_handle_t step = path_range.begin;
                while (true) {
                    auto next = graph.get_next_step(step);
                    dset.unite(id_to_entry[graph.get_id(graph.get_handle_of_step(step))],
                               id_to_entry[graph.get_id(graph.get_handle_of_step(next))]);
					if (next == path_range.end) break;
					step = next;
                }
            }
            // compute the block count and node to block mapping
            uint64_t n_dsets = 0;
            ska::flat_hash_map<uint64_t, uint64_t> dset_ids;
            ska::flat_hash_map<uint64_t, uint64_t> node_to_dset;
            for (auto& path_range : block.path_ranges) {
				step_handle_t step = path_range.begin;
                while (true) {
                    auto id_node = graph.get_id(graph.get_handle_of_step(step));
                    auto entry = id_to_entry[id_node];
                    auto id_dset = dset.find(entry);
                    auto f = dset_ids.find(id_dset);
                    if (f == dset_ids.end()) {
                        dset_ids[id_dset] = n_dsets++;
                    }
                    node_to_dset[id_node] = dset_ids[id_dset];
					if (step == path_range.end) break;
					step = graph.get_next_step(step);
                }
            }
            std::vector<block_t> blocks(n_dsets);
            // break the block apart
            for (auto& path_range : block.path_ranges) {
                step_handle_t step = path_range.begin;
                uint64_t dset_id = node_to_dset[graph.get_id(graph.get_handle_of_step(step))];
                while (true) {
                    uint64_t next_dset_id = node_to_dset[graph.get_id(graph.get_handle_of_step(step))];
                    assert(dset_id == next_dset_id);
					if (step == path_range.end) break;
					step = graph.get_next_step(step);
                }
                blocks[dset_id].path_ranges.push_back(path_range);
            }
            return blocks;
        };
    auto finalize_block =
        [&](block_t& block, std::vector<handle_t>& block_handles) {
            // collect the steps on all handles
            std::vector<step_handle_t> traversals;
            for (auto& handle : block_handles) {
                graph.for_each_step_on_handle(
                    handle,
                    [&](const step_handle_t& step) {
                        if (!seen_step(step)) {
                            traversals.push_back(step);
                        }
                    });
            }

			std::cout << "TRAVERSALS: " << std::endl;
			for (auto& traversal : traversals) {
				std::cout << "    node_id: " << graph.get_id(graph.get_handle_of_step(traversal)) << std::endl;
				std::cout << "    path_name: " << graph.get_path_name(graph.get_path_handle_of_step(traversal)) << std::endl;
				std::cout << "    path_position: " << step_index.get_position(traversal, graph) << std::endl;
			}

            std::vector<handle_t>().swap(block_handles);

            // sort them
            std::sort(
                traversals.begin(), traversals.end(),
                [&](const step_handle_t& a, const step_handle_t& b) {
                    return graph.get_path_name(graph.get_path_handle_of_step(a)) < graph.get_path_name(graph.get_path_handle_of_step(b)) || as_integer(graph.get_path_handle_of_step(a)) == as_integer(graph.get_path_handle_of_step(b)) && step_index.get_position(a, graph) < step_index.get_position(b, graph);
                });

			std::cout << "SORTED TRAVERSALS: " << std::endl;
			for (auto& traversal : traversals) {
				std::cout << "    node_id: " << graph.get_id(graph.get_handle_of_step(traversal)) << std::endl;
				std::cout << "    path_name: " << graph.get_path_name(graph.get_path_handle_of_step(traversal)) << std::endl;
				std::cout << "    path_position: " << step_index.get_position(traversal, graph) << std::endl;
			}
            // determine the path ranges in the block
            // break them when we pass some threshold for how much block-external sequence to include
            // (this parameter is meant to allow us to reduce dispersed collapses in the graph)
            // break them when they jump more than max_edge_jump in our graph sort order
            // TODO explore breaking when we have a significant change in coverage relative to the average in the region
            std::vector<path_range_t> path_ranges;
            for (auto& step : traversals) {
                if (path_ranges.empty()) {
                    path_ranges.push_back({step, step, 0});
                } else {
                    auto& path_range = path_ranges.back();
                    auto& last = path_range.end;
                    if (as_integer(graph.get_path_handle_of_step(last)) != as_integer(graph.get_path_handle_of_step(step))
                        || (step_index.get_position(step, graph)
                            - (step_index.get_position(last, graph) + graph.get_length(graph.get_handle_of_step(last)))
                            > max_path_jump)) {
                        // make a new range
                        path_ranges.push_back({step, step, 0});
                    } else {
                        // extend the range
                        last = step;
                    }
                }
            }
            // break the path ranges on seen steps
            for (auto& path_range : path_ranges) {
                uint64_t included_path_length = 0;
                // update the path range end to point to the one-past element
				// TODO FIX THIS
				if (graph.has_next_step(path_range.end)) {
					path_range.end = graph.get_next_step(path_range.end);
				}
                path_range_t* curr_path_range = nullptr;
                step_handle_t curr_step = path_range.begin;
                while (true) {
                    //if (curr_step == graph.path_end(graph.get_path_handle_of_step(curr_step))
                    if (curr_path_range == nullptr) {
                        block.path_ranges.emplace_back();
                        curr_path_range = &block.path_ranges.back();
                        curr_path_range->begin = curr_step;
                    }
                    curr_path_range->end = curr_step;
                    if (seen_step(curr_step)) {
                        curr_path_range = nullptr;
                    }
					if (curr_step == path_range.end) break;
					curr_step = graph.get_next_step(curr_step);
				}
                if (curr_path_range != nullptr) {
                    curr_path_range->end = curr_step;
                }
            }

            // erase any empty path ranges that we picked up
            // and any path ranges that are shorter than our minimum subpath size
            block.path_ranges.erase(
                std::remove_if(
                    block.path_ranges.begin(), block.path_ranges.end(),
                    [&graph](const path_range_t& path_range) {
                        return path_range.begin == path_range.end;
                    }),
                block.path_ranges.end());

            // finally, mark which steps we've kept and record the total length
            uint64_t _total_path_length = 0; // recalculate how much sequence we have
            //block.max_path_length = 0; // and the longest path range
            for (auto& path_range : block.path_ranges) {
                auto& included_path_length = path_range.length;
                included_path_length = 0;
                // here we need to break when we see significant nonlinearities
				step_handle_t curr_step = path_range.begin;
                while (true) {
                    mark_step(curr_step);
                    included_path_length += graph.get_length(graph.get_handle_of_step(curr_step));
					if (curr_step == path_range.end) break;
					curr_step = graph.get_next_step(curr_step);
				}
                _total_path_length += included_path_length;
            }

            if (_total_path_length > 0) {
                // order the path ranges from longest/shortest to shortest/longest
                // this gets called lots of times... probably best to make it std::sort or not parallel
                std::sort(
                        block.path_ranges.begin(), block.path_ranges.end(),
                        order_paths_from_longest
                        ?
                        [](const path_range_t& a,
                           const path_range_t& b) {
                            return a.length > b.length;
                        }
                        :
                        [](const path_range_t& a,
                           const path_range_t& b) {
                            return a.length < b.length;
                        }
                );
                // split blocks by graph topology
                // here weakly connected components of the graph are split apart
                // so that we do not compress disparate parts of the graph in one POA block
                for (auto& split : toposplit_block(block)) {
                    blockset.add_block(rank++, split);
                }
            };

            std::vector<path_range_t>().swap(block.path_ranges);
        };
    //uint64_t id = 0;
    std::stringstream blocks_banner;
    blocks_banner << "[smoothxg::smoothable_blocks] computing blocks for "
                    << graph.get_node_count() << " handles:";
    progress_meter::ProgressMeter blocks_progress(graph.get_node_count(), blocks_banner.str());

    uint64_t total_path_length = 0;
    ska::flat_hash_map<path_handle_t, std::pair<uint64_t, uint64_t>> path_coverage;

    graph.for_each_handle(
        [&](const handle_t& handle) {
            // how much sequence would we be adding to the block?
            int64_t handle_length = graph.get_length(handle);
            uint64_t sequence_to_add = 0;
            graph.for_each_step_on_handle(
                handle,
                [&](const step_handle_t& step) {
                    if (!seen_step(step)) {
                        sequence_to_add += handle_length;
                    }
                });
            // nb. this doesn't count duplicate traversals, but doing so would require running a full block traversal at every handle
            uint64_t max_path_length = 0;
            for (auto& p : path_coverage) {
                uint64_t path_len_est = std::round((double)p.second.first
                                                   / (p.second.second < block_handles.size()
                                                      ? 1.0
                                                      : (double) p.second.second / (double) block_handles.size()));
                max_path_length = std::max(path_len_est + handle_length, max_path_length);
            }
			/* FIXME
			std::cerr << "NODE_ID: " << graph.get_id(handle) << std::endl;
			std::cerr << "max_path_length: " << max_path_length << std::endl;
			*/
            // for each edge, find the jump length
            int64_t longest_edge_jump = 0;
            int64_t handle_vec_offset = node_offsets[graph.get_id(handle)];
			// FIXME
			// std::cerr << "handle_vec_offset: " << handle_vec_offset << std::endl;

			//int64_t handle_length = graph.get_length(handle);
            graph.follow_edges(
                handle, false,
                [&](const handle_t& o) {
                    int64_t other_vec_offset = node_offsets[graph.get_id(o)]
                        + (graph.get_is_reverse(o) ? graph.get_length(o) : 0);
                    int64_t jump = std::abs(other_vec_offset - (handle_vec_offset + handle_length));
                    longest_edge_jump = std::max(longest_edge_jump, jump);
                });
            graph.follow_edges(
                handle, true,
                [&](const handle_t& o) {
                    int64_t other_vec_offset = node_offsets[graph.get_id(o)]
                        + (graph.get_is_reverse(o) ? 0 : graph.get_length(o));
                    int64_t jump = std::abs(other_vec_offset - handle_vec_offset);
                    longest_edge_jump = std::max(longest_edge_jump, jump);
                });

            if (
                    // if it is not the first handle in the block
                    !block_handles.empty() &&

                    // if we add to the current block, do we go over our total path length?
                    // do we have a path that seems to exceed our length limit?
                    ((total_path_length + sequence_to_add > max_block_weight)
                     || (max_edge_jump && longest_edge_jump > max_edge_jump)
                     || max_path_length > max_block_path_length)
            ) {
                /*
                std::cerr << "block over weight "
                          << block.total_path_length << " " << sequence_to_add << " " << max_block_weight << std::endl;
                */

                // if so, finalize the last block and add the new one
                finalize_block(block, block_handles);

                total_path_length = sequence_to_add;
                path_coverage.clear();

            } else {
                // if not, add and update

                total_path_length += sequence_to_add;
            }

			// FIXME
			// std::cerr << "total_path_length: " << total_path_length << std::endl;

			graph.for_each_step_on_handle(
                handle,
                [&](const step_handle_t& step) {
                    if (!seen_step(step)) {
                        path_coverage[graph.get_path_handle_of_step(step)].first += handle_length;
                        path_coverage[graph.get_path_handle_of_step(step)].second++;
                    }
                });

            block_handles.push_back(handle);

            blocks_progress.increment(1);
        });

    blocks_progress.finish();

    if (block.path_ranges.empty()) {
        finalize_block(block, block_handles);
    }

    // at the end, we'll be left with some fragments of paths that aren't included in any blocks
    // that's ok, but we should see how much of a problem it is / should they be compressed?

    blockset.index(num_threads);
}

}
