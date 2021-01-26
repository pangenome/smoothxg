#include "utils.hpp"

namespace smoothxg {

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

void graph_deep_copy(odgi::graph_t* source,
                     odgi::graph_t* target) {
    // copy the now-compacted graph to our output_graph
    source->for_each_handle(
        [&](const handle_t& old_handle) {
            target->create_handle(
                source->get_sequence(old_handle),
                source->get_id(old_handle));
        });

    source->for_each_handle(
        [&](const handle_t& curr) {
            source->follow_edges(
                curr, false,
                [&](const handle_t& next) {
                    target->create_edge(
                        target->get_handle(source->get_id(curr),
                                                 source->get_is_reverse(curr)),
                        target->get_handle(source->get_id(next),
                                                 source->get_is_reverse(next)));
                });
            source->follow_edges(
                curr, true,
                [&](const handle_t& prev) {
                    target->create_edge(
                        target->get_handle(source->get_id(prev),
                                                 source->get_is_reverse(prev)),
                        target->get_handle(source->get_id(curr),
                                                 source->get_is_reverse(curr)));
                });
        });

    source->for_each_path_handle(
        [&](const path_handle_t& old_path) {
            path_handle_t new_path = target->create_path_handle(source->get_path_name(old_path));
            source->for_each_step_in_path(old_path, [&](const step_handle_t& step) {
                    handle_t old_handle = source->get_handle_of_step(step);
                    handle_t new_handle = target->get_handle(
                        source->get_id(old_handle),
                        source->get_is_reverse(old_handle));
                    target->append_step(new_path, new_handle);
                });
        });
}

}
