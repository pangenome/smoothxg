#include "prep.hpp"

namespace smoothxg {

// prep the graph into a given GFA file
// we'll then build the xg index on top of that in low memory

void prep_for_smoothing(
    const std::string& gfa_in,
    const std::string& gfa_out,
    const uint64_t& max_node_length
    ) {


    // load it into an odgi
    odgi::graph_t graph;
    odgi::gfa_to_handle(gfa_in, &graph);

    // chop it
    odgi::algorithms::chop(graph, max_node_length);

    // sort it using a short sorting pipeline
    // first toposort
    graph.apply_ordering(odgi::algorithms::topological_order(
                             &graph, true, false, false),
                         true);
    
    // then path-guided SGD

    // parameters that we might like to set
    uint64_t path_sgd_iter_max = 30; //args::get(p_sgd_iter_max) ? args::get(p_sgd_iter_max) : 30;
    double path_sgd_zipf_theta = 0.99; // args::get(p_sgd_zipf_theta) ? args::get(p_sgd_zipf_theta) : 0.99;
    double path_sgd_eps = 0.01; // args::get(p_sgd_eps) ? args::get(p_sgd_eps) : 0.01;
    double path_sgd_delta = 0; //args::get(p_sgd_delta) ? args::get(p_sgd_delta) : 0;
    std::set<std::string> path_sgd_use_paths;

    // path length interrogation
    std::function<uint64_t(const std::set<std::string> &,
                           const xp::XP &)> get_sum_path_lengths
            = [&](const std::set<std::string> &path_sgd_use_paths, const xp::XP &path_index) {
                uint64_t sum_path_length = 0;
                for (auto path_name : path_sgd_use_paths) {
                    handlegraph::path_handle_t path = path_index.get_path_handle(path_name);
                    sum_path_length += path_index.get_path_length(path);
                }
                return sum_path_length;
              };

    std::function<uint64_t(const std::set<std::string> &,
                           const xp::XP &)> get_max_path_length
            = [&](const std::set<std::string> &path_sgd_use_paths, const xp::XP &path_index) {
                uint64_t max_path_length = 0;
                for (auto path_name : path_sgd_use_paths) {
                    handlegraph::path_handle_t path = path_index.get_path_handle(path_name);
                    max_path_length = std::max(max_path_length, path_index.get_path_length(path));
                }
                return max_path_length;
              };

    xp::XP path_index;
    path_index.from_handle_graph(graph);

    uint64_t sum_path_length = get_sum_path_lengths(path_sgd_use_paths, path_index);
    float p_sgd_min_term_updates = 1; // -G parameter to odgi sort
    uint64_t path_sgd_min_term_updates = p_sgd_min_term_updates * sum_path_length;
    uint64_t path_sgd_zipf_space = get_max_path_length(path_sgd_use_paths, path_index);
    std::string path_sgd_seed = "pangenomic!";
    
    auto order
        = odgi::algorithms::path_linear_sgd_order(
            graph,
            path_index,
            path_sgd_use_paths,
            path_sgd_iter_max,
            path_sgd_min_term_updates,
            path_sgd_delta,
            path_sgd_eps,
            path_sgd_zipf_theta,
            path_sgd_zipf_space,
            odgi::get_thread_count(),
            false,
            path_sgd_seed);

    graph.apply_ordering(order, true);

    // groom
    odgi::graph_t groomed;
    odgi::algorithms::groom(graph, groomed);
    graph = groomed;

    // final toposort
    graph.apply_ordering(odgi::algorithms::topological_order(
                             &graph, true, false, false),
                         true);

    std::ofstream f(gfa_out);
    graph.to_gfa(f);

}

}
