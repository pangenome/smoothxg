#pragma once

#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <IITree.h> // cgranges
#include "xg.hpp"
//#include "index.hpp"

namespace smoothxg {

using namespace handlegraph;

// wrap some long-winded but simple functions for manipulating handles

inline uint64_t handle_rank(const handle_t& handle) {
    return number_bool_packing::unpack_number(handle);
}

inline bool handle_is_rev(const handle_t& handle) { 
    return number_bool_packing::unpack_bit(handle);
}

inline handle_t make_handle(const uint64_t& rank, const bool& is_rev) {
    return number_bool_packing::pack(rank, is_rev);
}

inline handle_t flip_handle(const handle_t& handle) {
    return number_bool_packing::pack(handle_rank(handle), !handle_is_rev(handle));
}

// hard coded ID/handle mapping
inline nid_t to_id(const handle_t& handle) {
    return handle_rank(handle) + 1;
}

inline handle_t max_handle(void) {
    return make_handle(std::numeric_limits<int64_t>::max(), false);
}

// seq_pos_t defines an oriented position in the sequence space of the graph
// 
// the position can be on the forward or reverse complement of the graph
// 
// to simplify use during clustering, the offset is always measured
// relative to the beginning of the sequence vector on that strand
// 
// the encoding (with orientation in the most significant bit)
// allows us to increment and decrement positions using standard operators
//
typedef std::uint64_t seq_pos_t;

// these are helper functions for creating and interrogating the seq_pos_t's
struct seq_pos {
    constexpr static uint64_t OFFSET_BITS      = 63;
    constexpr static uint64_t ORIENTATION_MASK = static_cast<uint64_t>(1) << OFFSET_BITS;
    constexpr static uint64_t OFFSET_MASK      = ORIENTATION_MASK - 1;
    static seq_pos_t encode(const uint64_t& offset, const bool& reverse_complement) {
        return offset | (reverse_complement ? ORIENTATION_MASK : 0);
    }
    static bool is_rev(const seq_pos_t& pos) { return pos & ORIENTATION_MASK; }
    static uint64_t offset(const seq_pos_t& pos) { return pos & OFFSET_MASK; }
    static std::string to_string(const seq_pos_t& pos) {
        std::stringstream ss;
        ss << offset(pos) << (is_rev(pos) ? "-" : "+");
        return ss.str();
    }
};

struct anchor_t {
    seq_pos_t query_begin = 0;
    seq_pos_t query_end = 0;
    seq_pos_t target_begin = 0;
    seq_pos_t target_end = 0;
    double max_chain_score = 0;
    anchor_t* best_predecessor = nullptr;
    anchor_t(const seq_pos_t& qb,
             const seq_pos_t& qe,
             const seq_pos_t& tb,
             const seq_pos_t& te)
        : query_begin(qb)
        , query_end(qe)
        , target_begin(tb)
        , target_end(te) { }
};

/*
std::vector<anchor_t> anchors_for_query(const gyeet_index_t& index,
                                        const char* seq,
                                        const size_t& len);
*/

struct chain_t {
    std::vector<anchor_t*> anchors;
    double score = 0;
    double mapping_quality = std::numeric_limits<double>::min();
    bool is_secondary = false;
    bool processed(void) {
        return mapping_quality != std::numeric_limits<double>::min();
    }
    // inner target boundaries
    seq_pos_t target_begin = 0;
    seq_pos_t target_end = 0;
    // query boundaries are fixed
    seq_pos_t query_begin(void) const { return anchors.front()->query_begin; }
    seq_pos_t query_end(void) const { return anchors.back()->query_end; }
    void compute_boundaries(const uint64_t& seed_length, const double& mismatch_rate);
};

struct chain_node_t {
    chain_t* chain = nullptr;
    chain_node_t* best_predecessor = nullptr;
    double max_superchain_score = 0;
    bool used = false;
    chain_node_t(chain_t* c) : chain(c) { }
};

std::vector<chain_t>
chains(std::vector<anchor_t>& anchors,
       const uint64_t& seed_length,
       const uint64_t& max_gap,
       const double& mismatch_rate,
       const uint64_t& chain_min_n_anchors,
       const uint64_t bandwidth = 50,
       const double secondary_chain_threshold = 0.5,
       const double max_mapq = 60);

double score_anchors(const anchor_t& a,
                     const anchor_t& b,
                     const uint64_t& seed_length,
                     const uint64_t& max_gap);

uint64_t chain_query_length(const chain_t& chain);

struct superchain_t {
    std::vector<chain_t*> chains;
    double score = 0;
    bool is_secondary = true;
    /*
    double mapping_quality = std::numeric_limits<double>::min();
    //bool is_secondary = false;
    bool processed(void) {
        return mapping_quality != std::numeric_limits<double>::min();
    }
    */
};

std::vector<superchain_t>
superchains(std::vector<chain_t>& chains,
            const uint64_t& kmer_length,
            const double& mismatch_rate,
            const double& chain_overlap_max,
            const uint64_t bandwidth = 1000);

double score_chain_nodes(const chain_node_t& a,
                         const chain_node_t& b,
                         const uint64_t& kmer_length,
                         const double& overlap_max);

std::vector<std::vector<superchain_t>> collinear_blocks(const xg::XG& graph);

}