#include "breaks.hpp"

namespace smoothxg {

using namespace handlegraph;

// break the path ranges at likely VNTR boundaries
// and break the path ranges to be shorter than our "max" sequence size input to spoa
void break_blocks(const xg::XG& graph,
                  std::vector<block_t>& blocks,
                  const uint64_t& max_poa_length,
                  const uint64_t& min_copy_length,
                  const uint64_t& max_copy_length,
                  const uint64_t& min_autocorr_z,
                  const uint64_t& autocorr_stride,
                  const bool& order_paths_from_longest) {

    const VectorizableHandleGraph& vec_graph = dynamic_cast<const VectorizableHandleGraph&>(graph);

    std::cerr << "[smoothxg::break_blocks] cutting blocks that contain sequences longer than max-poa-length (" << max_poa_length << ")" << std::endl;

    uint64_t n_cut_blocks = 0;
    uint64_t n_repeat_blocks = 0;
    for (auto& block : blocks) {
        // check if we have sequences that are too long
        bool to_break = false;
        for (auto& path_range : block.path_ranges) {
            if (path_range.length > max_poa_length) {
                to_break = true;
                break;
            }
        }
        if (!to_break) continue; // skip if we're spoa-able
        // otherwise let's see if we've got repeats that we can use to chop things up
        // find if there is a repeat
        std::vector<sautocorr::repeat_t> repeats;
        for (auto& path_range : block.path_ranges) {
            // steps in id space
            std::string seq;
            std::string name = graph.get_path_name(graph.get_path_handle_of_step(path_range.begin));
            for (step_handle_t step = path_range.begin;
                 step != path_range.end;
                 step = graph.get_next_step(step)) {
                seq.append(graph.get_sequence(graph.get_handle_of_step(step)));
            }
            if (seq.length() < 2*min_copy_length) continue;
            //std::cerr << "on " << name << "\t" << seq.length() << std::endl;
            std::vector<uint8_t> vec(seq.begin(), seq.end());
            sautocorr::repeat_t result = sautocorr::repeat(vec,
                                                           min_copy_length,
                                                           max_copy_length,
                                                           min_copy_length,
                                                           min_autocorr_z,
                                                           autocorr_stride);
            repeats.push_back(result);
            /*
            std::cerr << name
                      << "\t" << seq.length()
                      << "\t" << result.length
                      << "\t" << result.z_score << std::endl;
            */
        }
        // if there is, set the cut length to some fraction of it
        std::vector<double> lengths;
        for (auto& repeat : repeats) {
            if (repeat.length > 0) {
                lengths.push_back(repeat.length);
            }
        }
        uint64_t cut_length;
        bool found_repeat = !lengths.empty();
        if (found_repeat) {
            double repeat_length = sautocorr::vec_mean(lengths.begin(), lengths.end());
            cut_length = std::round(repeat_length / 2.0);
            ++n_repeat_blocks;
            //std::cerr << "found repeat of " << repeat_length << " cutting to " << cut_length << std::endl;
        } else {
            // if not, chop blindly
            cut_length = max_poa_length;
        }
        ++n_cut_blocks;
        std::vector<path_range_t> chopped_ranges;
        for (auto& path_range : block.path_ranges) {

            if (!found_repeat && path_range.length < cut_length) {
                chopped_ranges.push_back(path_range);
                continue;
            }
            // now find outlier clusters based on stdev and mean
            // extract a minimum viable repeat length
            // scan across the step vector, looking for where the repeat region begins and ends
            // cut at the repeat boundaries

            // Q: should we determine the repeat length for each sequence or all?
            // each is simple, but maybe expensive
            // all could provide higher precision, but it's muddier

            // if this doesn't work, we're going to blindly cut anyway
            uint64_t last_cut = 0;
            step_handle_t last_end = path_range.begin;
            //path_range_t* new_range = nullptr;
            uint64_t pos = 0;
            step_handle_t step;
            for (step = path_range.begin;
                 step != path_range.end;
                 step = graph.get_next_step(step)) {
                //handle_t h = graph.get_handle_of_step(step);
                //uint64_t id = graph.get_id(h);
                //int64_t node_pos = vec_graph.node_vector_offset(id);
                pos += graph.get_length(graph.get_handle_of_step(step));
                if (pos - last_cut > cut_length) {
                    step_handle_t next = graph.get_next_step(step);
                    chopped_ranges.push_back({last_end, next, pos - last_cut});
                    last_end = next;
                    last_cut = pos;
                }
            }
            if (step != last_end) {
                chopped_ranges.push_back({last_end, step, pos - last_cut});
            }
        }
        block.path_ranges = chopped_ranges;
        // order the path ranges from longest/shortest to shortest/longest
        ips4o::parallel::sort(
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
        block.broken = true;
        block.is_repeat = found_repeat;
    }
    std::cerr << "[smoothxg::break_blocks] cut " << n_cut_blocks << " blocks of which " << n_repeat_blocks << " had repeats" << std::endl;
}

}
