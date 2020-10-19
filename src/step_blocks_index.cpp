#include "step_blocks_index.hpp"

namespace smoothxg {

ska::flat_hash_map<std::string, IITree<uint64_t , uint64_t>> generate_step_rank_to_path_ranges_trees(
    std::vector<smoothxg::block_t>& blocks) {
    ska::flat_hash_map<std::string, IITree<uint64_t , uint64_t>> happy_trees_map;
    // TODO PLAY AROUND WITH THE MAP
    for (uint64_t i; i < blocks.size(); i++) {
        auto& block = blocks[i];
        // find out each path's step_handle ranges of a block
        // record this in the happy tress map
        for (auto path_range : block.path_ranges) {
            path_range_t* path_range_pointer = &path_range;
            // TODO go from pointer to integer
            uint64_t magic_pointer;
            // TODO go from step_handle_t to path
            std::string path_name;
            // TODO find out if there already is an entry for *path_name* in the map
            // TODO if so, fetch the entry tree, else create a new one
            IITree<uint64_t , uint64_t> happy_tree;
            // TODO go from step_handle_t to step_rank
            uint64_t begin_rank;
            uint64_t end_rank;
            happy_tree.add(begin_rank, end_rank, magic_pointer);
        }
    }
    return happy_trees_map;
}

}
