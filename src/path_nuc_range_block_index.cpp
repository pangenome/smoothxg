#include "path_nuc_range_block_index.hpp"

namespace smoothxg {

ska::flat_hash_map<std::string, IITree<uint64_t , uint64_t>> generate_path_nuc_range_block_index(
    std::vector<smoothxg::block_t>& blocks, xg::XG& graph) {
    ska::flat_hash_map<std::string, IITree<uint64_t , uint64_t>> happy_trees_map;
    /*
    // TODO PLAY AROUND WITH THE MAP
    IITree<uint64_t, uint64_t> some_tree = happy_trees_map["a"];
    if (some_tree.size() == 0) {
        std::cerr << "SO IT GOES" << std::endl;
    }
    */
    for (uint64_t i; i < blocks.size(); i++) {
        auto& block = blocks[i];
        // find out each path's step_handle ranges of a block
        // record this in the happy tress map
        for (auto path_range : block.path_ranges) {
            // go from step_handle_t to path name
            std::string path_name = graph.get_path_name(graph.get_path_handle_of_step(path_range.begin));
            // find out if there already is an entry for *path_name* in the map
            // if so, fetch the entry tree, else create a new one
            IITree<uint64_t , uint64_t> happy_tree = happy_trees_map[path_name];
            if (happy_tree.size() == 0) {
                IITree<uint64_t , uint64_t> happy_tree_;
                happy_tree_.add(path_range.nuc_begin, path_range.nuc_end + 1, i);
                happy_trees_map[path_name] = happy_tree_;
            } else {
                happy_tree.add(path_range.nuc_begin, path_range.nuc_end + 1, i);
            }
            /*
            std::cerr << path_range.nuc_begin << std::endl;
            std::cerr << path_range.nuc_end << std::endl;
            */
        }
    }
    // std::cerr << "happy_trees_map.size() " << happy_trees_map.size() << std::endl;
    for (auto& happy_tree : happy_trees_map) {
        happy_tree.second.index();
        /*
        std::cerr << "happy_tree.size() " << happy_tree.second.size() << std::endl;
        std::vector<size_t> a;
        IITree tree = happy_tree.second;
        tree.overlap(19, 119, a); // retrieve overlaps
        for (size_t i = 0; i < a. size(); i++) {
            std::cerr << tree.start(a[i]) << std::endl;
            std::cerr << tree.end(a[i]) << std::endl;
            std::cerr << tree.data(a[i]) << std::endl;
        }
        */
    }
    /*
            IITree<int, int> tree;
tree.add(12, 34, 0); // add an interval
tree.add(0, 23, 1);
tree.add(34, 56, 2);
tree.index(); // index
std::vector<size_t> a;
tree.overlap(23, 25, a); // retrieve overlaps
for (size_t i = 0; i < a.size(); ++i) {
    printf("%d\t%d\t%d\n", tree.start(a[i]), tree.end(a[i]), tree.data(a[i]));
    std::cerr << tree.start(a[i]) << std::endl;
    std::cerr << tree.end(a[i]) << std::endl;
    std::cerr << tree.data(a[i]) << std::endl;
}
     */

    return happy_trees_map;
}

}
