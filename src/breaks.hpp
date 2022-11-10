#pragma once

#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <vector>
#include <vector>
#include <numeric>
#include <cmath>
#include "WFA/edit/edit_cigar.hpp"
#include "WFA/gap_affine/affine_wavefront_align.hpp"
#include "WFA/utils/commons.hpp"
#include "blocks.hpp"
#include "sautocorr.hpp"
#include "xg.hpp"
#include "odgi/dna.hpp"

namespace smoothxg {

using namespace handlegraph;

// break the path ranges at likely VNTR boundaries
// and break the path ranges to be shorter than our "max" sequence size input to spoa
void break_blocks(const xg::XG& graph,
                  blockset_t*& blockset,
                  const double &length_ratio_min,
                  const uint64_t& min_length_mash_based_clustering,
                  const double& block_group_identity,
                  const double& block_group_est_identity,
                  const uint64_t& kmer_size,
                  const uint64_t& min_dedup_depth_for_block_splitting,
                  const uint64_t& min_dedup_depth_for_mash_clustering,
                  const uint64_t& max_poa_length,
                  const uint64_t& min_copy_length,
                  const uint64_t& max_copy_length,
                  const uint64_t& min_autocorr_z,
                  const uint64_t& autocorr_stride,
                  const bool& order_paths_from_longest,
                  const bool& break_repeats,
                  const uint64_t& thread_count,
                  const bool& write_block_to_split_fastas,
				  const std::string& smoothxg_iter);

}
