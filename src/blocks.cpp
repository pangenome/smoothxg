#include "blocks.hpp"

namespace smoothxg {

std::vector<block_t>
smoothable_blocks(
    const xg::XG& graph,
    const uint64_t& max_block_weight,
    const uint64_t& max_path_jump) {
    // iterate over the handles in their vectorized order
    std::vector<block_t> blocks;
    std::vector<std::vector<bool>> seen_steps;
    graph.for_each_path_handle(
        [&](const path_handle_t& path) {
            seen_steps.emplace_back();
            seen_steps.back().resize(graph.get_step_count(path));
        });
    auto seen_step =
        [&](const step_handle_t& step) {
            // in xg, the first half of the step is the path handle, which is it's rank + 1
            // and the second half of the step is the rank in the path
            return seen_steps[as_integers(step)[0]-1][as_integers(step)[1]];
        };
    auto mark_step =
        [&](const step_handle_t& step) {
            seen_steps[as_integers(step)[0]-1][as_integers(step)[1]] = true;
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
            for (auto& step : traversals) {
                if (block.path_ranges.empty()) {
                    block.path_ranges.push_back({step, step});
                } else {
                    auto& path_range = block.path_ranges.back();
                    auto& last = path_range.end;
                    if (path_rank(last) != path_rank(step)
                        || (graph.get_position_of_step(step)
                            - (graph.get_position_of_step(last) + graph.get_length(graph.get_handle_of_step(last)))
                            > max_path_jump)) {
                        // make a new range
                        block.path_ranges.push_back({step, step});
                    }
                }
            }
            // mark which steps we've kept
            for (auto& path_range : block.path_ranges) {
                step_handle_t curr = path_range.begin;
                mark_step(curr);
                while (curr != path_range.end) {
                    curr = graph.get_next_step(curr);
                    mark_step(curr);
                }
            }
        };
    graph.for_each_handle(
        [&](const handle_t& handle) {
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
                    // if so, finalize the last block and add the new one
                    finalize_block(block);
                    blocks.emplace_back();
                    blocks.back().handles.push_back(handle);
                } else {
                    // if not, add and update
                    block.handles.push_back(handle);
                }
            }
        });
    if (blocks.back().path_ranges.empty()) {
        finalize_block(blocks.back());
    }
    // at the end, we'll be left with some fragments of paths that aren't included in any blocks
    // that's ok, but we should see how much of a problem it is / should they be compressed?
    return blocks;
}

}
