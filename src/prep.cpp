#include "prep.hpp"

namespace smoothxg {

// prep the graph into a given GFA file
// we'll then build the xg index on top of that in low memory

void prep(
    const std::string& gfa_in,
    const std::string& gfa_out,
    const uint64_t& max_node_length,
    const float& p_sgd_min_term_updates) {

    // load it into an odgi
    odgi::graph_t graph;
    odgi::gfa_to_handle(gfa_in, &graph);

    // chop it
    odgi::algorithms::chop(graph, max_node_length);

    // sort it using a short sorting pipeline
    // first toposort
    graph.apply_ordering(odgi::algorithms::topological_order(
                             &graph, true, false, true),
                         true);
    
    // then path-guided SGD

    // parameters that we might like to set
    uint64_t path_sgd_iter_max = 30; //args::get(p_sgd_iter_max) ? args::get(p_sgd_iter_max) : 30;
    double path_sgd_zipf_theta = 0.99; // args::get(p_sgd_zipf_theta) ? args::get(p_sgd_zipf_theta) : 0.99;
    double path_sgd_eps = 0.01; // args::get(p_sgd_eps) ? args::get(p_sgd_eps) : 0.01;
    double path_sgd_delta = 0; //args::get(p_sgd_delta) ? args::get(p_sgd_delta) : 0;
    std::vector<handlegraph::path_handle_t> path_sgd_use_paths;
    graph.for_each_path_handle(
        [&](const handlegraph::path_handle_t &path) {
            path_sgd_use_paths.push_back(path);
        });

    // path length interrogation
    auto get_sum_path_lengths
        = [&](const std::vector<handlegraph::path_handle_t> &path_sgd_use_paths, const xp::XP &path_index) {
              uint64_t sum_path_length = 0;
              for (auto& path : path_sgd_use_paths) {
                  sum_path_length += path_index.get_path_length(path);
              }
              return sum_path_length;
          };
    auto get_max_path_length
        = [&](const std::vector<handlegraph::path_handle_t> &path_sgd_use_paths, const xp::XP &path_index) {
              uint64_t max_path_length = 0;
              for (auto& path : path_sgd_use_paths) {
                  max_path_length = std::max(max_path_length, path_index.get_path_length(path));
              }
              return max_path_length;
          };

    xp::XP path_index;
    path_index.from_handle_graph(graph);

    uint64_t sum_path_length = get_sum_path_lengths(path_sgd_use_paths, path_index);
    uint64_t path_sgd_min_term_updates = p_sgd_min_term_updates * graph.get_node_count();
    uint64_t path_sgd_zipf_space = get_max_path_length(path_sgd_use_paths, path_index);
    double path_sgd_max_eta = graph.get_node_count();
    std::string path_sgd_seed = "pangenomic!";

    /*
    std::cerr
        << path_sgd_iter_max << " "
        << path_sgd_min_term_updates << " "
        << path_sgd_delta << " "
        << path_sgd_eps << " "
        << path_sgd_zipf_theta << " "
        << path_sgd_zipf_space << std::endl;
    */
/*
    std::vector<handle_t> path_linear_sgd_order(const PathHandleGraph &graph,
                                            const xp::XP &path_index,
                                            const std::vector<path_handle_t>& path_sgd_use_paths,
                                            const uint64_t &iter_max,
                                            const uint64_t &iter_with_max_learning_rate,
                                            const uint64_t &min_term_updates,
                                            const double &delta,
                                            const double &eps,
                                            const double &eta_max,
                                            const double &theta,
                                            const uint64_t &space,
                                            const uint64_t &nthreads,
                                            const bool &progress,
                                            const std::string &seed,
                                            const bool &snapshot,
                                            std::vector<std::vector<handle_t>> &snapshots);
*/

    std::vector<std::vector<handlegraph::handle_t>> null_snapshots;
    
    auto order
        = odgi::algorithms::path_linear_sgd_order(
            graph,
            path_index,
            path_sgd_use_paths,
            path_sgd_iter_max,
            0,
            path_sgd_min_term_updates,
            path_sgd_delta,
            path_sgd_eps,
            path_sgd_max_eta,
            path_sgd_zipf_theta,
            path_sgd_zipf_space,
            odgi::get_thread_count(),
            true,
            path_sgd_seed,
            false,
            null_snapshots);

    graph.apply_ordering(order, true);

    // groom
    odgi::graph_t groomed;
    odgi::algorithms::groom(graph, groomed);
    graph = groomed;

    // final toposort
    graph.apply_ordering(odgi::algorithms::topological_order(
                             &graph, true, false, true),
                         true);

    std::ofstream f(gfa_out);
    graph.to_gfa(f);

}

}
