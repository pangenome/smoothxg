#include "prep.hpp"
#include <filesystem>

namespace smoothxg {

#define MAX_NUMBER_OF_ZIPF_DISTRIBUTIONS 100

// prep the graph into a given GFA file
// we'll then build the xg index on top of that in low memory

void prep(
    const std::string& gfa_in,
    const std::string& gfa_out,
    const uint64_t& max_node_length,
    const float& p_sgd_min_term_updates,
    const bool& toposort,
    const std::string& basename,
    const uint64_t& num_threads,
	const std::string& smoothxg_iter) {

    // load it into an odgi
    odgi::graph_t graph;

    std::filesystem::path graph_path = gfa_in;
    std::cerr << smoothxg_iter << "::prep] Loading graph for prep " << gfa_in << std::endl;

    if (graph_path.extension() == ".gfa") {
        odgi::gfa_to_handle(gfa_in, &graph, true, num_threads, true);
    } else {
        ifstream f(gfa_in.c_str());
        graph.deserialize(f);
        f.close();
    }
    graph.set_number_of_threads(num_threads);

    // sort it using a short sorting pipeline equivalent to `odgi sort -p Ygs`

    // first path-guided SGD

    // parameters that we might like to set
    uint64_t path_sgd_iter_max = 100; //args::get(p_sgd_iter_max) ? args::get(p_sgd_iter_max) : 30;
    double path_sgd_zipf_theta = 0.99; // args::get(p_sgd_zipf_theta) ? args::get(p_sgd_zipf_theta) : 0.99;
    double path_sgd_eps = 0.01; // args::get(p_sgd_eps) ? args::get(p_sgd_eps) : 0.01;
    double path_sgd_delta = 0; //args::get(p_sgd_delta) ? args::get(p_sgd_delta) : 0;
    double path_sgd_cooling = 0.5; // initiate cooling halfway through our iterations
    std::vector<handlegraph::path_handle_t> path_sgd_use_paths;
    graph.for_each_path_handle(
        [&](const handlegraph::path_handle_t &path) {
            path_sgd_use_paths.push_back(path);
        });

    // path length interrogation
    std::function<uint64_t(const std::vector<handlegraph::path_handle_t> &,
                           const xp::XP &)> get_sum_path_step_count
        = [&](const std::vector<handlegraph::path_handle_t> &path_sgd_use_paths, const xp::XP &path_index) {
              uint64_t sum_path_step_count = 0;
              for (auto& path : path_sgd_use_paths) {
                  sum_path_step_count += path_index.get_path_step_count(path);
              }
              return sum_path_step_count;
          };
    std::function<uint64_t(const std::vector<handlegraph::path_handle_t> &,
                           const xp::XP &)> get_max_path_step_count
        = [&](const std::vector<handlegraph::path_handle_t> &path_sgd_use_paths, const xp::XP &path_index) {
              uint64_t max_path_step_count = 0;
              for (auto& path : path_sgd_use_paths) {
                  max_path_step_count = std::max(max_path_step_count, path_index.get_path_step_count(path));
              }
              return max_path_step_count;
          };

	std::function<uint64_t(const std::vector<handlegraph::path_handle_t> &,
						   const xp::XP &)> get_max_path_length
			= [&](const std::vector<handlegraph::path_handle_t> &path_sgd_use_paths, const xp::XP &path_index) {
				uint64_t max_path_length = std::numeric_limits<uint64_t>::min();
				for (auto &path : path_sgd_use_paths) {
					max_path_length = std::max(max_path_length, path_index.get_path_length(path));
				}
				return max_path_length;
			};

    std::cerr << smoothxg_iter << "::prep] building path index" << std::endl;
    xp::XP path_index;
    path_index.from_handle_graph(graph, basename, num_threads);

    uint64_t sum_path_step_count = get_sum_path_step_count(path_sgd_use_paths, path_index);
    uint64_t path_sgd_min_term_updates = p_sgd_min_term_updates * sum_path_step_count;
    uint64_t max_path_step_count = get_max_path_step_count(path_sgd_use_paths, path_index);
    uint64_t path_sgd_zipf_space = get_max_path_length(path_sgd_use_paths, path_index); //std::min((uint64_t)10000, max_path_step_count);
    double path_sgd_max_eta = max_path_step_count * max_path_step_count;
    uint64_t path_sgd_zipf_space_max = 100;
	uint64_t path_sgd_zipf_max_number_of_distributions = std::max((uint64_t) path_sgd_zipf_space_max + 1,
																  (uint64_t) MAX_NUMBER_OF_ZIPF_DISTRIBUTIONS);
	std::cerr << smoothxg_iter << "::prep] path_sgd_zipf_space_max: " << path_sgd_zipf_space_max << std::endl;
	std::cerr << smoothxg_iter << "::prep] path_sgd_zipf_max_number_of_distributions: " << path_sgd_zipf_max_number_of_distributions << std::endl;
    uint64_t path_sgd_zipf_space_quantization_step = 100;
	if (path_sgd_zipf_space > path_sgd_zipf_space_max && path_sgd_zipf_max_number_of_distributions > path_sgd_zipf_space_max) {
		path_sgd_zipf_space_quantization_step = std::max(
				(uint64_t) 2,
				(uint64_t) ceil( (double) (path_sgd_zipf_space - path_sgd_zipf_space_max) / (double) (path_sgd_zipf_max_number_of_distributions - path_sgd_zipf_space_max))
		);
	}

    std::string path_sgd_seed = "pangenomic!";

    uint64_t path_sgd_iter_max_learning_rate = 0; // don't use this max iter stuff
    const std::string snapshot_prefix;
	const bool target_sorting = false;
	std::vector<bool> target_nodes;

    std::cerr << smoothxg_iter << "::prep] sorting graph" << std::endl;
    auto order
        = odgi::algorithms::path_linear_sgd_order(
            graph,
            path_index,
            path_sgd_use_paths,
            path_sgd_iter_max,
            path_sgd_iter_max_learning_rate,
            path_sgd_min_term_updates,
            path_sgd_delta,
            path_sgd_eps,
            path_sgd_max_eta,
            path_sgd_zipf_theta,
            path_sgd_zipf_space,
            path_sgd_zipf_space_max,
            path_sgd_zipf_space_quantization_step,
            path_sgd_cooling,
            num_threads,
            true,
            path_sgd_seed,
            false,
            snapshot_prefix,
            false,
            snapshot_prefix,
			target_sorting,
			target_nodes);

    graph.set_number_of_threads(num_threads);
    graph.apply_ordering(order, true);

    // groom
	const std::vector<handlegraph::path_handle_t> target_paths;
    graph.apply_ordering(odgi::algorithms::groom(graph, true, target_paths));
    graph.set_number_of_threads(num_threads);

    // final toposort
    if (toposort) {
        graph.apply_ordering(odgi::algorithms::topological_order(
                                 &graph, true, true, true),
                             true);
    }

    std::cerr << smoothxg_iter << "::prep] chopping graph to " << max_node_length << std::endl;
    // chop it (preserves order)
    odgi::algorithms::chop(graph, max_node_length, num_threads, true);

    std::cerr << smoothxg_iter << "::prep] writing graph " << gfa_out << std::endl;

    graph_path = gfa_out;
    std::ofstream f(gfa_out);
    if (graph_path.extension() == ".gfa") {
        graph.to_gfa(f);
    } else {
        graph.serialize(f);
    }
    f.close();
}

}
