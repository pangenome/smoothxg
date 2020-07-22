#include "blocks.hpp"

namespace smoothxg {

std::vector<block_t>
smoothable_blocks(
    const xg::XG& graph,
    const uint64_t& max_block_weight,
    const uint64_t& max_path_jump,
    const uint64_t& min_subpath) {
    // iterate over the handles in their vectorized order
    std::vector<block_t> blocks;
    std::vector<std::vector<bool>> seen_steps;
    std::cerr << "[smoothxg::smoothable_blocks] computing blocks" << std::endl;
    graph.for_each_path_handle(
        [&](const path_handle_t& path) {
            seen_steps.emplace_back();
            seen_steps.back().resize(graph.get_step_count(path));
        });
    auto seen_step =
        [&](const step_handle_t& step) {
            // in xg, the first half of the step is the path handle, which is it's rank + 1
            // and the second half of the step is the rank in the path
            return seen_steps[path_rank(step)-1][step_rank(step)];
        };
    auto mark_step =
        [&](const step_handle_t& step) {
            seen_steps[path_rank(step)-1][step_rank(step)] = true;
        };
    auto finalize_block =
        [&](block_t& block) {
            // collect the steps on all handles
            std::vector<step_handle_t> traversals;
            for (auto& handle : block.handles) {
                graph.for_each_step_on_handle(
                    handle,
                    [&](const step_handle_t& step) {
                        if (!seen_step(step)) {
                            traversals.push_back(step);
                        }
                    });
            }
            // sort them
            ips4o::parallel::sort(
                traversals.begin(), traversals.end(),
                [&](const step_handle_t& a, const step_handle_t& b) {
                    return path_rank(a) < path_rank(b) || path_rank(a) == path_rank(b) && step_rank(a) < step_rank(b);
                });
            // determine the path ranges in the block
            // break them when we pass some threshold for how much block-external sequence to include
            // (this parameter is meant to allow us to reduce dispersed collapses in the graph)
            std::vector<path_range_t> path_ranges;
            for (auto& step : traversals) {
                if (path_ranges.empty()) {
                    path_ranges.push_back({step, step});
                } else {
                    auto& path_range = path_ranges.back();
                    auto& last = path_range.end;
                    if (path_rank(last) != path_rank(step)
                        || (graph.get_position_of_step(step)
                            - (graph.get_position_of_step(last) + graph.get_length(graph.get_handle_of_step(last)))
                            > max_path_jump)) {
                        // make a new range
                        path_ranges.push_back({step, step});
                    } else {
                        // extend the range
                        last = step;
                    }
                }
            }
            // break the blocks on ranges of seen steps
            //block.path_ranges
            for (auto& path_range : path_ranges) {
                uint64_t included_path_length = 0;
                // update the path range end to point to the one-past element
                path_range.end = graph.get_next_step(path_range.end);
                path_range_t* curr_path_range = nullptr;
                step_handle_t curr_step;
                for (curr_step = path_range.begin;
                     curr_step != path_range.end;
                     curr_step = graph.get_next_step(curr_step)) { 
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
                    [&graph,&min_subpath](const path_range_t& path_range) {
                        uint64_t range_length =
                            graph.get_position_of_step(graph.get_previous_step(path_range.end))
                            - graph.get_position_of_step(path_range.begin);
                        return (min_subpath ? range_length < min_subpath : false)
                            || path_range.begin == path_range.end;
                    }),
                block.path_ranges.end());

            // finally, mark which steps we've kept and record the total length
            block.total_path_length = 0; // recalculate how much sequence we have
            block.max_path_length = 0; // and the longest path range
            for (auto& path_range : block.path_ranges) {
                auto& included_path_length = path_range.length;
                included_path_length = 0;
                /*
                std::cerr << "on path range for " << graph.get_path_name(graph.get_path_handle_of_step(path_range.begin))
                          << " " << graph.get_id(graph.get_handle_of_step(path_range.begin))
                          << "-"
                          << graph.get_id(graph.get_handle_of_step(graph.get_previous_step(path_range.end))) << std::endl;
                */
                for (step_handle_t curr_step = path_range.begin;
                     curr_step != path_range.end;
                     curr_step = graph.get_next_step(curr_step)) {
                    //std::cerr << "on step " << graph.get_id(graph.get_handle_of_step(curr_step)) << std::endl;
                    mark_step(curr_step);
                    included_path_length += graph.get_length(graph.get_handle_of_step(curr_step));
                }
                block.total_path_length += included_path_length;
                block.max_path_length = std::max(included_path_length,
                                                 block.max_path_length);
            }

            // order the path ranges from longest to shortest
            ips4o::parallel::sort(
                block.path_ranges.begin(), block.path_ranges.end(),
                [](const path_range_t& a,
                   const path_range_t& b) {
                    return a.length > b.length;
                });
            /*
            std::cerr << "block----" << std::endl;
            for (auto& path_range : block.path_ranges) {
                std::cerr << "path_range " << path_range.length << " "
                          << graph.get_path_name(graph.get_path_handle_of_step(path_range.begin))
                          << " " << graph.get_id(graph.get_handle_of_step(path_range.begin))
                          << "-"
                          << graph.get_id(graph.get_handle_of_step(graph.get_previous_step(path_range.end))) << std::endl;
            }
            */                    
        };
    //uint64_t id = 0;
    graph.for_each_handle(
        [&](const handle_t& handle) {
            if (graph.get_id(handle) % 100 == 0) {
                std::cerr << std::fixed << std::showpoint << std::setprecision(3)
                          << "[smoothxg::smoothable_blocks] computing blocks "
                          << (float)graph.get_id(handle) / (float)graph.get_node_count() * 100 << "%\r";
            }
            if (blocks.empty()) {
                blocks.emplace_back();
                auto& block = blocks.back();
                block.handles.push_back(handle);
            } else {
                // how much sequence would we be adding to the block?
                uint64_t handle_length = graph.get_length(handle);
                uint64_t sequence_to_add = 0;
                graph.for_each_step_on_handle(
                    handle,
                    [&](const step_handle_t& step) {
                        if (!seen_step(step)) {
                            sequence_to_add += handle_length;
                        }
                    });
                auto& block = blocks.back();
                // if we add to the current block, do we go over our total path length?
                if (block.total_path_length + sequence_to_add > max_block_weight) {
                    /*
                    std::cerr << "block over weight "
                              << block.total_path_length << " " << sequence_to_add << " " << max_block_weight << std::endl;
                    */
                    // if so, finalize the last block and add the new one
                    finalize_block(block);
                    blocks.emplace_back();
                    blocks.back().handles.push_back(handle);
                } else {
                    // if not, add and update
                    block.handles.push_back(handle);
                    block.total_path_length += sequence_to_add;
                }
            }
        });
    std::cerr << "[smoothxg::smoothable_blocks] computing blocks 100.00%" << std::endl;
    if (blocks.back().path_ranges.empty()) {
        finalize_block(blocks.back());
    }
    // at the end, we'll be left with some fragments of paths that aren't included in any blocks
    // that's ok, but we should see how much of a problem it is / should they be compressed?
    return blocks;
}

}
