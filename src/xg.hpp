#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <queue>
#include <omp.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <functional>
#include <utility>
#include <tuple>
#include <sys/types.h>
#include <dirent.h>

#include "sdsl/bit_vectors.hpp"
#include "sdsl/enc_vector.hpp"
#include "sdsl/dac_vector.hpp"
#include "sdsl/vlc_vector.hpp"
#include "sdsl/wavelet_trees.hpp"
#include "sdsl/csa_wt.hpp"
#include "sdsl/suffix_arrays.hpp"

#include "handlegraph/types.hpp"
#include "handlegraph/iteratee.hpp"
#include "handlegraph/util.hpp"
#include "handlegraph/handle_graph.hpp"
#include "handlegraph/path_position_handle_graph.hpp"

#include "mmmultimap.hpp"


namespace xg {

using namespace handlegraph;
using nid_t = handlegraph::nid_t;

typedef std::tuple<nid_t, bool, size_t> pos_t;
inline pos_t make_pos_t(const nid_t& id, bool is_rev, const size_t& off) {
    return std::make_tuple(id, is_rev, off);
}

/// Extract the id of the node a pos_t is on.
inline nid_t id(const pos_t& pos) {
    return std::get<0>(pos);
}

/// Return true if a pos_t is on the reverse strand of its node.
inline bool is_rev(const pos_t& pos) {
    return std::get<1>(pos);
}

/// Get the offset along the selected strand of the node from a pos_t.
inline size_t offset(const pos_t& pos) {
    return std::get<2>(pos);
}

/// Get a reference to the Node ID of a pos_t.
inline nid_t& get_id(pos_t& pos) {
    return std::get<0>(pos);
}

/// Get a reference to the reverse flag of a pos_t.
inline bool& get_is_rev(pos_t& pos) {
    return std::get<1>(pos);
}

/// Get a reference to the offset field of a pos_t, which counts along the selected strand of the node.
inline size_t& get_offset(pos_t& pos) {
    return std::get<2>(pos);
}

/// Return true if a pos_t is unset.
inline bool is_empty(const pos_t& pos) {
    return (id(pos) == 0);
}

/// Reverse a pos_t and get a pos_t at the same **point between bases**, going the other direction.
/// To get a pos_t to the same *base*, subtract 1 from the resulting offset or call reverse_base_pos().
inline pos_t reverse(const pos_t& pos, size_t node_length) {
    pos_t rev = pos;
    // swap the offset onto the other strand
    get_offset(rev) = node_length - offset(rev);
    // invert the position
    get_is_rev(rev) = !is_rev(rev);
    return rev;
}

/// Reverse a pos_t and get a pos_t at the same **base**, going the other direction.
inline pos_t reverse_base_pos(const pos_t& pos, size_t node_length) {
    pos_t rev = pos;
    // swap the offset onto the other strand
    get_offset(rev) = (node_length - 1) - offset(rev);
    // invert the position
    get_is_rev(rev) = !is_rev(rev);
    return rev;
}

/// Print a pos_t to a stream.
inline std::ostream& operator<<(std::ostream& out, const pos_t& pos) {
    return out << id(pos) << (is_rev(pos) ? "-" : "+") << offset(pos);
}

static const char complement[256] = {'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 8
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 16
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 24
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 32
                                     'N', 'N', 'N', '$', '#', 'N', 'N', 'N', // 40 GCSA stop/start characters
                                     'N', 'N', 'N', 'N', 'N', '-', 'N', 'N', // 48
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 56
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 64
                                     'N', 'T', 'V', 'G', 'H', 'N', 'N', 'C', // 72
                                     'D', 'N', 'N', 'M', 'N', 'K', 'N', 'N', // 80
                                     'N', 'Q', 'Y', 'W', 'A', 'A', 'B', 'S', // 88
                                     'N', 'R', 'N', 'N', 'N', 'N', 'N', 'N', // 96
                                     'N', 't', 'v', 'g', 'h', 'N', 'N', 'c', // 104
                                     'd', 'N', 'N', 'm', 'N', 'k', 'n', 'N', // 112
                                     'N', 'q', 'y', 'w', 'a', 'a', 'b', 's', // 120
                                     'N', 'r', 'N', 'N', 'N', 'N', 'N', 'N', // 128
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 136
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 144
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 152
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 160
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 168
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 176
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 184
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 192
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 200
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 208
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 216
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 224
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 232
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 240
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 248
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};// 256

inline char reverse_complement(const char& c) {
    return complement[c];
}

inline std::string reverse_complement(const std::string& seq) {
    std::string rc;
    rc.assign(seq.rbegin(), seq.rend());
    for (auto& c : rc) {
        c = complement[c];
    }
    return rc;
}
    
inline void reverse_complement_in_place(std::string& seq) {
    size_t swap_size = seq.size() / 2;
    for (size_t i = 0, j = seq.size() - 1; i < swap_size; i++, j--) {
        char tmp = seq[i];
        seq[i] = complement[seq[j]];
        seq[j] = complement[tmp];
    }
    
    if (seq.size() % 2) {
        seq[swap_size] = complement[seq[swap_size]];
    }
}


class XGPath;

/**
 * Thrown when attempting to interpret invalid data as an XG index.
 */
class XGFormatError : public std::runtime_error {
    // Use the runtime_error constructor
    using std::runtime_error::runtime_error;
};

/**
 * Provides succinct storage for a graph, its positional paths, and a set of
 * embedded threads.
 */
class XG : public PathPositionHandleGraph, public SerializableHandleGraph, public VectorizableHandleGraph {
public:
    
    ////////////////////////////////////////////////////////////////////////////
    // Here are the ways we can construct XG objects (from graph data or files)
    ////////////////////////////////////////////////////////////////////////////
    
    XG(void) = default;
    ~XG(void);

    // We cannot move, assign, or copy until we add code to point SDSL suppots
    // at the new addresses for their vectors.
    XG(const XG& other) = delete;
    XG(XG&& other) = delete;
    XG& operator=(const XG& other) = delete;
    XG& operator=(XG&& other) = delete;

    /// build the graph from another simple graph
    void from_handle_graph(const HandleGraph& graph);

    /// build the graph from another path handle graph
    void from_path_handle_graph(const PathHandleGraph& graph);

    /// Use external enumerators to drive graph construction
    void from_enumerators(const std::function<void(const std::function<void(const std::string& seq, const nid_t& node_id)>&)>& for_each_sequence,
                          const std::function<void(const std::function<void(const nid_t& from, const bool& from_rev,
                                                                            const nid_t& to, const bool& to_rev)>&)>& for_each_edge,
                          const std::function<void(const std::function<void(const std::string& path_name,
                                                                            const nid_t& node_id, const bool& is_rev,
                                                                            const std::string& cigar, const bool& is_empty,
                                                                            const bool& is_circular)>&)>& for_each_path_element,
                          bool validate = false, std::string basename = "");

    /// Use a memory-mapped GFA file to build the index in low memory
    void from_gfa(const std::string& gfa_filename, bool validate = false, std::string basename = "");

    void to_gfa(std::ostream& out) const;

    // Use an existing handle graph to build the index
    //void from_path_handle_graph(const PathHandleGraph& other);
    
    // What version of an XG is this designed to read?
    const static uint32_t CURRENT_VERSION = 13;
               
    // Load this XG index from a stream. Throw an XGFormatError if the stream
    // does not produce a valid XG file.
    void load(std::istream& in);
    
    // Alias for load() to match the SerializableHandleGraph interface
    void deserialize(std::istream& in);
    
    void serialize(std::ostream& out) const;
    size_t serialize_and_measure(std::ostream& out, sdsl::structure_tree_node* s = nullptr, std::string name = "") const;
    
    ////////////////////////////////////////////////////////////////////////////
    // Basic API
    ////////////////////////////////////////////////////////////////////////////
    
    // General public statisitcs
    size_t seq_length = 0;
    size_t node_count = 0;
    size_t edge_count = 0;
    size_t path_count = 0;
    
    ////////////////////////////////////////////////////////////////////////////
    // Here is the handle graph API
    ////////////////////////////////////////////////////////////////////////////
    
    /// Look up the handle for the node with the given ID in the given orientation
    virtual handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const;
    /// Get the ID from a handle
    virtual nid_t get_id(const handle_t& handle) const;
    /// If the node with the given id exists
    virtual bool has_node(nid_t node_id) const;
    /// Get the orientation of a handle
    virtual bool get_is_reverse(const handle_t& handle) const;
    /// Invert the orientation of a handle (potentially without getting its ID)
    virtual handle_t flip(const handle_t& handle) const;
    /// Get the length of a node
    virtual size_t get_length(const handle_t& handle) const;
    /// Get the sequence of a node, presented in the handle's local forward
    /// orientation.
    virtual std::string get_sequence(const handle_t& handle) const;
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue.
    virtual bool follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const;
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee returns false.
    virtual bool for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    /// Return the number of nodes in the graph
    virtual size_t get_node_count() const;
    /// Get the minimum node ID used in the graph, if any are used
    virtual nid_t min_node_id() const;
    /// Get the maximum node ID used in the graph, if any are used
    virtual nid_t max_node_id() const;
    /// Returns one base of a handle's sequence, in the orientation of the
    /// handle.
    virtual char get_base(const handle_t& handle, size_t index) const;
    /// Returns a substring of a handle's sequence, in the orientation of the
    /// handle. If the indicated substring would extend beyond the end of the
    /// handle's sequence, the return value is truncated to the sequence's end.
    virtual std::string get_subsequence(const handle_t& handle, size_t index, size_t size) const;
    
    // TODO: There's currently no really good efficient way to implement
    // get_degree; we have to decode each edge to work out what node side it is
    // on. So we use the default implementation.
    
    ////////////////////////
    // Path handle graph API
    ////////////////////////
   
    /// Returns the number of paths stored in the graph
    size_t get_path_count(void) const;
    /// Determine if a path with a given name exists
    bool has_path(const std::string& path_name) const;
    /// Look up the path handle for the given path name
    path_handle_t get_path_handle(const std::string& path_name) const;
    /// Look up the name of a path from a handle to it
    std::string get_path_name(const path_handle_t& path_handle) const;
    /// Look up whether a path is circular
    bool get_is_circular(const path_handle_t& path_handle) const;
    /// Returns the number of node steps in the path
    size_t get_step_count(const path_handle_t& path_handle) const;
    /// Returns the total length of sequence in the path
    size_t get_path_length(const path_handle_t& path_handle) const;
    /// Get a node handle (node ID and orientation) from a handle to a step on a path
    handle_t get_handle_of_step(const step_handle_t& step_handle) const;
    /// Returns a handle to the path that an step is on
    path_handle_t get_path_handle_of_step(const step_handle_t& step_handle) const;
    /// Get a handle to the first step, or in a circular path to an arbitrary step
    /// considered "first". If the path is empty, returns the past-the-last step
    /// returned by path_end.
    step_handle_t path_begin(const path_handle_t& path_handle) const;
    /// Get a handle to a fictitious position past the end of a path. This position is
    /// return by get_next_step for the final step in a path in a non-circular path.
    /// Note that get_next_step will *NEVER* return this value for a circular path.
    step_handle_t path_end(const path_handle_t& path_handle) const;
    /// Get a handle to the last step, which will be an arbitrary step in a circular path that
    /// we consider "last" based on our construction of the path. If the path is empty
    /// then the implementation must return the same value as path_front_end().
    step_handle_t path_back(const path_handle_t& path_handle) const;
    /// Get a handle to a fictitious position before the beginning of a path. This position is
    /// return by get_previous_step for the first step in a path in a non-circular path.
    /// Note: get_previous_step will *NEVER* return this value for a circular path.
    step_handle_t path_front_end(const path_handle_t& path_handle) const;
    /// Returns true if the step is not the last step in a non-circular path.
    bool has_next_step(const step_handle_t& step_handle) const;
    /// Returns true if the step is not the first step in a non-circular path.
    bool has_previous_step(const step_handle_t& step_handle) const;
    /// Returns a handle to the next step on the path. If the given step is the final step
    /// of a non-circular path, returns the past-the-last step that is also returned by
    /// path_end. In a circular path, the "last" step will loop around to the "first" (i.e.
    /// the one returned by path_begin).
    /// Note: to iterate over each step one time, even in a circular path, consider
    /// for_each_step_in_path.
    step_handle_t get_next_step(const step_handle_t& step_handle) const;
    /// Returns a handle to the previous step on the path. If the given step is the first
    /// step of a non-circular path, this method has undefined behavior. In a circular path,
    /// it will loop around from the "first" step (i.e. the one returned by path_begin) to
    /// the "last" step.
    /// Note: to iterate over each step one time, even in a circular path, consider
    /// for_each_step_in_path.
    step_handle_t get_previous_step(const step_handle_t& step_handle) const;
    
    /// Execute a function on each path in the graph
    bool for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const;
    /// Executes a function on each step of a handle in any path.
    bool for_each_step_on_handle_impl(const handle_t& handle, const std::function<bool(const step_handle_t&)>& iteratee) const;
    /// Unpack the path position and orientation information alongside the steps
    bool for_each_step_position_on_handle(const handle_t& handle, const std::function<bool(const step_handle_t&, const bool&, const uint64_t&)>& iteratee) const;
    /// Gets the position of a given step in the path it's from
    size_t get_position_of_step(const step_handle_t& step) const;
    /// Get the step at a given position
    step_handle_t get_step_at_position(const path_handle_t& path, const size_t& position) const;
    
    ////////////////////////////////////////////////////////////////////////////
    // Higher-level graph API
    ////////////////////////////////////////////////////////////////////////////
    /// id to rank helper function

    size_t id_to_rank(const nid_t& id) const;
    size_t handle_rank(const handle_t& handle) const;
    /// rank to id helper function
    nid_t rank_to_id(const size_t& rank) const;
    size_t max_node_rank(void) const;
    int64_t node_at_seq_pos(const size_t& pos) const;
    size_t node_vector_offset(const nid_t& id) const;
    nid_t node_at_vector_offset(const size_t& offset) const;
    size_t max_path_rank(void) const;
    size_t node_graph_idx(const nid_t& id) const;
    const XGPath& get_path(const std::string& name) const;
    std::string path_name(const size_t& rank) const;
    std::vector<size_t> path_ranks_by_prefix(const std::string& prefix) const;
    std::vector<std::string> path_names_by_prefix(const std::string& prefix) const;
    std::vector<path_handle_t> paths_of_handle(const handle_t& handle) const;
    bool path_contains_handle(const std::string& name, const handle_t& handle) const;
    std::pair<pos_t, int64_t> next_path_position(const pos_t& pos, const int64_t& max_search) const;
    std::pair<nid_t, std::vector<path_handle_t> > nearest_path_node(const nid_t& id, int64_t max_steps = 16) const;
    int64_t min_approx_path_distance(const nid_t& id1, const nid_t& id2) const;
    std::vector<size_t> position_in_path(const handle_t& handle, const path_handle_t& path) const;
    std::unordered_map<path_handle_t, std::vector<size_t> > position_in_paths(const handle_t& handle, const size_t& offset) const;
    void for_path_range(const std::string& name, int64_t start, int64_t stop,
                        std::function<void(const handle_t&)> lambda, bool is_rev) const;
    std::unordered_map<path_handle_t, std::vector<std::pair<size_t, bool> > > offsets_in_paths(const pos_t& gpos) const;
    std::unordered_map<path_handle_t, std::vector<std::pair<size_t, bool> > > nearest_offsets_in_paths(const pos_t& pos, int64_t max_search) const;
    handle_t handle_at_path_position(const path_handle_t& path, size_t pos) const;
    size_t node_start_at_path_position(const path_handle_t& path, size_t pos) const;
    pos_t graph_pos_at_path_position(const path_handle_t& path, size_t path_pos) const;
    pos_t graph_pos_at_path_position(const std::string& name, size_t path_pos) const;
    char pos_char(nid_t id, bool is_rev, size_t off) const;
    std::string pos_substr(nid_t id, bool is_rev, size_t off, size_t len) const;
    edge_t edge_from_encoding(const nid_t& from, const nid_t& to, int type) const;
    size_t edge_index(const edge_t& edge) const;
    size_t get_g_iv_size(void) const;

    char start_marker = '#';
    char end_marker = '$';

private:

    ////////////////////////////////////////////////////////////////////////////
    // Here is the New Way (locally traversable graph storage)
    // Everything should be rewritten in terms of these members
    ////////////////////////////////////////////////////////////////////////////

    /// locally traversable graph storage
    /// 
    /// Encoding designed for efficient compression, cache locality, and relativistic traversal of the graph.
    ///
    /// node := { header, edges_to, edges_from }
    /// header := { node_id, node_start, node_length, edges_to_count, edges_from_count }
    /// node_id := integer
    /// node_start := integer (offset in s_iv)
    /// node_length := integer
    /// edges_to_count := integer
    /// edges_from_count := integer
    /// edges_to := { edge_to, ... }
    /// edges_from := { edge_from, ... }
    /// edge_to := { offset_to_previous_node, edge_type }
    /// edge_to := { offset_to_next_node, edge_type }
    sdsl::int_vector<> g_iv;
    /// delimit node records to allow lookup of nodes in g_civ by rank
    sdsl::bit_vector g_bv;
    sdsl::rank_support_v<1> g_bv_rank;
    sdsl::bit_vector::select_1_type g_bv_select;
    
    // Let's define some offset ints
    const static int G_NODE_ID_OFFSET = 0;
    const static int G_NODE_SEQ_START_OFFSET = 1;
    const static int G_NODE_LENGTH_OFFSET = 2;
    const static int G_NODE_TO_COUNT_OFFSET = 3;
    const static int G_NODE_FROM_COUNT_OFFSET = 4;
    const static int G_NODE_HEADER_LENGTH = 5;

    const static int G_EDGE_OFFSET_OFFSET = 0;
    const static int G_EDGE_TYPE_OFFSET = 1;
    const static int G_EDGE_LENGTH = 2;

    int edge_type(const handle_t& from, const handle_t& to) const;
    
    /// This is a utility function for the edge exploration. It says whether we
    /// want to visit an edge depending on its type, whether we're the to or
    /// from node, whether we want to look left or right, and whether we're
    /// forward or reverse on the node.
    bool edge_filter(int type, bool is_to, bool want_left, bool is_reverse) const;
    
    // This loops over the given number of edge records for the given g node,
    // starting at the given start g vector position. For all the edges that are
    // wanted by edge_filter given the is_to, want_left, and is_reverse flags,
    // the iteratee is called. Returns true if the iteratee never returns false,
    // or false (and stops iteration) as soon as the iteratee returns false.
    bool do_edges(const size_t& g, const size_t& start, const size_t& count,
                  bool is_to, bool want_left, bool is_reverse, const std::function<bool(const handle_t&)>& iteratee) const;
    
    // Use memmapped indexing to construct the node-to-path indexes once
    // XGPath's have been created (used during construction)
    void index_node_to_path(const std::string& basename);
    
    ////////////////////////////////////////////////////////////////////////////
    // Here are the bits we need to keep around to talk about the sequence
    ////////////////////////////////////////////////////////////////////////////
    
    // sequence/integer vector
    sdsl::int_vector<> s_iv;
    // node starts in sequence, provides id schema
    // rank_1(i) = id
    // select_1(id) = i
    sdsl::bit_vector s_bv; // node positions in siv
    sdsl::rank_support_v<1> s_bv_rank;
    sdsl::bit_vector::select_1_type s_bv_select;
    
    ////////////////////////////////////////////////////////////////////////////
    // And here are the bits for tracking actual node IDs
    ////////////////////////////////////////////////////////////////////////////

    // maintain old ids from input, ranked as in s_iv and s_bv
    int64_t min_id = 0; // id ranges don't have to start at 0
    int64_t max_id = 0;
    sdsl::int_vector<> r_iv; // ids-id_min is the rank

    ////////////////////////////////////////////////////////////////////////////
    // Here is path storage
    ////////////////////////////////////////////////////////////////////////////

    // path names
    sdsl::int_vector<> pn_iv; // path names
    sdsl::csa_wt<> pn_csa; // path name compressed suffix array
    sdsl::bit_vector pn_bv;  // path name starts in uncompressed version of csa
    sdsl::rank_support_v<1> pn_bv_rank;
    sdsl::bit_vector::select_1_type pn_bv_select;
    sdsl::int_vector<> pi_iv; // path ids by rank in the path names

    std::vector<XGPath*> paths; // path structure

    // node->path membership
    sdsl::int_vector<> np_iv;
    // node->path rank
    sdsl::int_vector<> nr_iv;
    // node->path position/orientation
    sdsl::int_vector<> nx_iv;
    sdsl::bit_vector np_bv; // entity delimiters in both vectors
    //sdsl::rank_support_v<1> np_bv_rank;
    sdsl::bit_vector::select_1_type np_bv_select;
};

class XGPath {
public:
    XGPath(void) = default;
    ~XGPath(void) = default;
    // Path name is required here only for complaining intelligently when
    // something goes wrong. We can also spit out the total unique members,
    // because in here is the most efficient place to count them.
    XGPath(const std::string& path_name,
           const std::vector<handle_t>& path,
           bool is_circular,
           XG& graph);
    // Path names are stored in the XG object, in a compressed fashion, and are
    // not duplicated here.
    
    // These contain rank and select supports and so cannot move or be copied
    // without code to update them.
    XGPath(const XGPath& other) = delete;
    XGPath(XGPath&& other) = delete;
    XGPath& operator=(const XGPath& other) = delete;
    XGPath& operator=(XGPath&& other) = delete;
    handle_t min_handle;
    sdsl::enc_vector<> handles;
    //sdsl::rrr_vector directions; // forward or backward through nodes
    sdsl::rrr_vector<> offsets;
    sdsl::rrr_vector<>::rank_1_type offsets_rank;
    sdsl::rrr_vector<>::select_1_type offsets_select;
    bool is_circular = false;
    void load(std::istream& in);
    void load_from_old_version(std::istream& in, uint32_t file_version, const XG& graph);
    size_t serialize(std::ostream& out,
                     sdsl::structure_tree_node* v = NULL,
                     std::string name = "") const;

    // Get the node orientation at a 0-based offset.
    nid_t node(size_t offset) const;
    bool is_reverse(size_t offset) const;
    handle_t local_handle(const handle_t& handle) const;
    handle_t external_handle(const handle_t& handle) const;
    handle_t handle(size_t offset) const;
    handle_t handle_at_position(size_t pos) const;
    size_t handle_start(size_t offset) const;
    size_t step_rank_at_position(size_t pos) const;
};

/**
 * Temporary files. Create with create() and remove with remove(). All
 * temporary files will be deleted when the program exits normally or with
 * std::exit(). The files will be created in a directory determined from
 * environment variables, though this can be overridden with set_dir().
 * The interface is thread-safe.
 */
namespace temp_file {

/// Create a temporary file starting with the given base name
std::string create(const std::string& base);

/// Create a temporary file
std::string create();

/// Remove a temporary file
void remove(const std::string& filename);

/// Set a temp dir, overriding system defaults and environment variables.
void set_dir(const std::string& new_temp_dir);

/// Get the current temp dir
std::string get_dir();

} // namespace temp_file

}
