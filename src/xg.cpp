#include "xg.hpp"

//#include "ips4o.hpp"
#include "mmmultimap.hpp"

#include <bitset>
#include <arpa/inet.h>
#include <mutex>

#include <handlegraph/util.hpp>

#include "gfakluge.hpp"

//#define VERBOSE_DEBUG
//#define debug_algorithms
//#define debug_component_index
#define debug_path_index

namespace xg {

int dna3bit(char c) {
    switch (c) {
    case 'A':
        return 0;
    case 'T':
        return 1;
    case 'C':
        return 2;
    case 'G':
        return 3;
    default:
        return 4;
    }
}

char revdna3bit(int i) {
    switch (i) {
    case 0:
        return 'A';
    case 1:
        return 'T';
    case 2:
        return 'C';
    case 3:
        return 'G';
    default:
        return 'N';
    }
}

XG::~XG(void) {
    // Clean up any created XGPaths
    while (!paths.empty()) {
        delete paths.back();
        paths.pop_back();
    }
}
    
void XG::deserialize(std::istream& in) {
    // simple alias to match an external interface
    load(in);
}

void XG::load(std::istream& in) {

    if (!in.good()) {
        throw XGFormatError("Index file does not exist or index stream cannot be read");
    }

    // Version 0 is the last XG format without an explicit version specifier.
    // If we find a version specifier we will up this.
    uint32_t file_version = 0;

    // We need to look for the magic value
    char buffer;
    in.get(buffer);
    if (buffer == 'X') {
        in.get(buffer);
        if (buffer == 'G') {
            // We found the magic value!
            
            // Don't put it back, but the next 4 bytes are a version number.
            in.read((char*) &file_version, sizeof(file_version));
            // Make sure to convert from network to host byte order
            file_version = ntohl(file_version);
            
        } else {
            // Put back both characters
            in.unget();
            in.unget();
        }        
    } else {
        // Put back the one character
        in.unget();
    }
    
    if (file_version > CURRENT_VERSION) {
        // This XG file is from the distant future.
        throw XGFormatError("XG index file version " + std::to_string(file_version) +
                            " is too new to interpret (max version = " + std::to_string(CURRENT_VERSION) + ")");
    }
    
    try {
        
        ////////////////////////////////////////////////////////////////////////
        // DO NOT CHANGE THIS CODE without creating a new XG version:
        // 1. Increment OUTPUT_VERSION to a new integer.
        // 2. Change the serialization code.
        // 3. Add a new case here (or enhance an existing case) with new deserialization code.
        ////////////////////////////////////////////////////////////////////////
        switch (file_version) {
        
        case 5: // Fall through
        case 6:
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
            std::cerr << "warning:[XG] Loading an out-of-date XG format."
                      << "For better performance over repeated loads, consider recreating this XG index." << std::endl;
            // Fall through
        case 13:
            {
                sdsl::read_member(seq_length, in);
                sdsl::read_member(node_count, in);
                sdsl::read_member(edge_count, in);
                sdsl::read_member(path_count, in);
                size_t entity_count = node_count + edge_count;
                //cerr << sequence_length << ", " << node_count << ", " << edge_count << endl;
                sdsl::read_member(min_id, in);
                sdsl::read_member(max_id, in);
                
                if (file_version <= 8) {
                    // Load the old id int vector to skip
                    sdsl::int_vector<> i_iv;
                    i_iv.load(in);
                }
                r_iv.load(in);
                
                g_iv.load(in);
                g_bv.load(in);
                g_bv_rank.load(in, &g_bv);
                g_bv_select.load(in, &g_bv);

                s_iv.load(in);
                s_bv.load(in);
                s_bv_rank.load(in, &s_bv);
                s_bv_select.load(in, &s_bv);

                if (file_version <= 11) {
                    // Skip over gPBWT thread names
                    {
                        sdsl::csa_bitcompressed<> tn_csa;
                        tn_csa.load(in);
                    }
                    {
                        sdsl::sd_vector<> tn_cbv;
                        tn_cbv.load(in);
                        sdsl::sd_vector<>::rank_1_type tn_cbv_rank;
                        tn_cbv_rank.load(in, &tn_cbv);
                        sdsl::sd_vector<>::select_1_type tn_cbv_select;
                        tn_cbv_select.load(in, &tn_cbv);
                    }
                }
                
                if (file_version >= 7 && file_version <= 11) {
                    // There is a single haplotype count here for all components
                    // We ignore this part of the gPBWT
                    size_t haplotype_count;
                    sdsl::read_member(haplotype_count, in);
                }
                
                if (file_version <= 11) {
                    // Skip thread positions in gPBWT
                    {
                        sdsl::vlc_vector<> tin_civ;
                        tin_civ.load(in);
                    }
                    {
                        sdsl::vlc_vector<> tio_civ;
                        tio_civ.load(in);
                    }
                    {
                        sdsl::wt_int<> side_thread_wt;
                        side_thread_wt.load(in);
                    }
                }
                
                pn_iv.load(in);
                pn_csa.load(in);
                pn_bv.load(in);
                pn_bv_rank.load(in, &pn_bv);
                pn_bv_select.load(in, &pn_bv);
                pi_iv.load(in);
                sdsl::read_member(path_count, in);
                for (size_t i = 0; i < path_count; ++i) {
                    auto path = new XGPath;
                    // Load the path, giving it the file version and a
                    // rank-to-ID comversion function for format upgrade
                    // purposes.
                    if (file_version <= 12) {
                        path->load_from_old_version(in, file_version, *this);
                    }
                    else {
                        path->load(in);
                    }
                    paths.push_back(path);
                }
                
                if (file_version <= 12) {
                    // skip over the old node-to-path indexes
                    {
                        sdsl::int_vector<> old_np_iv;
                        old_np_iv.load(in);
                    }
                    {
                        sdsl::bit_vector old_np_bv;
                        old_np_bv.load(in);
                        sdsl::rank_support_v<1> old_np_bv_rank;
                        old_np_bv_rank.load(in, &np_bv);
                        sdsl::bit_vector::select_1_type old_np_bv_select;
                        old_np_bv_select.load(in, &np_bv);
                    }
                    
                    // create the new node-to-path indexes
                    index_node_to_path(temp_file::create());
                }
                else {
                    // we're in the more recent encoding, so we can load
                    // the node-to-path indexes directly
                    
                    np_bv.load(in);
                    np_bv_select.load(in, &np_bv);
                    np_iv.load(in);
                    nr_iv.load(in);
                    nx_iv.load(in);
                }
                
                if (file_version >= 6 && file_version <= 10) {
                    // load and ignore the component path set indexes (which have
                    // now been exported)
                    {
                        sdsl::int_vector<> path_ranks_iv;
                        path_ranks_iv.load(in);
                    }
                    {
                        sdsl::bit_vector path_ranks_bv;
                        path_ranks_bv.load(in);
                    }
                }
                
                if (file_version <= 11) {
                    // load and ignore the gPBWT entity vectors
                    {
                        sdsl::vlc_vector<> h_civ;
                        h_civ.load(in);
                    }
                    // and the thread starts
                    {
                        sdsl::vlc_vector<> ts_civ;
                        ts_civ.load(in);
                    }
                    // and the B arrays
                    {
                        sdsl::wt_rlmn<sdsl::sd_vector<>> bs_single_array;
                        bs_single_array.load(in);
                    }
                }
            }
            break;
        default:
            throw XGFormatError("Unimplemented XG format version: " + std::to_string(file_version));
        }
    } catch (const XGFormatError& e) {
        // Pass XGFormatErrors through
        throw e;
    } catch (const std::bad_alloc& e) {
        // We get std::bad_alloc generally if we try to read arbitrary data as an xg index.
        throw XGFormatError("XG input data not in XG version " + std::to_string(file_version) + " format (" + e.what() + ")");
    } catch (const std::exception& e) {
        // Other things will get re-thrown with a hint.
        std::cerr << "error [xg]: Unexpected error parsing XG data. Is it in version " << file_version << " XG format?" << std::endl;
        throw e;
    }
}

void XGPath::load(std::istream& in) {
    sdsl::read_member(min_handle, in);
    handles.load(in);
    offsets.load(in);
    offsets_rank.load(in, &offsets);
    offsets_select.load(in, &offsets);
    sdsl::read_member(is_circular, in);    
}

void XGPath::load_from_old_version(std::istream& in, uint32_t file_version, const XG& graph) {
    
    if (file_version < 8) {
        // skip over some members from early versions
        {
            sdsl::rrr_vector<> nodes;
            nodes.load(in);
        }
        {
            sdsl::rrr_vector<>::rank_1_type nodes_rank;
            nodes_rank.load(in);
        }
        {
            sdsl::rrr_vector<>::select_1_type nodes_select;
            nodes_select.load(in);
        }
    }
    
    // convert the old ID and direction vectors into the handles vector
    {
        // first make an int vector of handles
        sdsl::int_vector<> handles_iv;
        {
            int64_t id_offset = 0;
            if (file_version >= 8) {
                // IDs are stored relative to a minimum offset
                sdsl::read_member(id_offset, in);
            }
            
            sdsl::wt_gmr<> ids;
            ids.load(in);
            
            sdsl::sd_vector<> directions;
            directions.load(in);
            // compute the minimum handle
            min_handle = handlegraph::as_handle(numeric_limits<int64_t>::max());
            for (size_t i = 0; i < ids.size(); ++i) {
                min_handle = as_handle(min(as_integer(min_handle),
                                           as_integer(graph.get_handle(ids[i] + id_offset - 1, directions[i]))));
            }
            // convert the vector into a handle-based one with a min handle offset
            sdsl::util::assign(handles_iv, sdsl::int_vector<>(ids.size()));
            for (size_t i = 0; i < ids.size(); ++i) {
                handles_iv[i] = as_integer(graph.get_handle(ids[i] + id_offset - 1, directions[i])) - as_integer(min_handle);
            }
        }
        
        // re-encode the handles int vector with a variable-length encoding
        sdsl::util::assign(handles, sdsl::enc_vector<>(handles_iv));
    }
    
    // skip the rank and position vectors
    {
        sdsl::int_vector<> ranks;
        ranks.load(in);
    }
    {
        sdsl::int_vector<> positions;
        positions.load(in);
    }
    
    // the offsets vectors have the same encoding
    offsets.load(in);
    offsets_rank.load(in, &offsets);
    offsets_select.load(in, &offsets);
    
    if (file_version >= 10) {
        // there is support for circular paths in this version
        sdsl::read_member(is_circular, in);
    }
    else {
        // previous versions are interpreted as not circular
        is_circular = false;
    }
}

size_t XGPath::serialize(std::ostream& out,
                         sdsl::structure_tree_node* v,
                         std::string name) const {
    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_t written = 0;
    written += sdsl::write_member(min_handle, out, child, "min_handle" + name);
    written += handles.serialize(out, child, "path_handles_" + name);
    written += offsets.serialize(out, child, "path_node_starts_" + name);
    written += offsets_rank.serialize(out, child, "path_node_starts_rank_" + name);
    written += offsets_select.serialize(out, child, "path_node_starts_select_" + name);
    written += sdsl::write_member(is_circular, out, child, "is_circular_" + name);
    sdsl::structure_tree::add_size(child, written);
    return written;
}

XGPath::XGPath(const std::string& path_name,
               const std::vector<handle_t>& path,
               bool is_circular,
               XG& graph) {

    // The circularity flag is just a normal bool
    this->is_circular = is_circular;

    // node ids, the literal path
    sdsl::int_vector<> handles_iv;
    sdsl::util::assign(handles_iv, sdsl::int_vector<>(path.size()));
    // directions of traversal (typically forward, but we allow backwards
    sdsl::bit_vector directions_bv;
    sdsl::util::assign(directions_bv, sdsl::bit_vector(path.size()));

    size_t path_off = 0;
    size_t members_off = 0;
    size_t positions_off = 0;
    size_t path_length = 0;
    min_handle = as_handle(std::numeric_limits<int64_t>::max());

    // determine min node id
    for (size_t i = 0; i < path.size(); ++i) {
        min_handle = as_handle(std::min(as_integer(min_handle), as_integer(path[i])));
    }
    // determine total length and record node ids
    for (size_t i = 0; i < path.size(); ++i) {
        const handle_t& handle = path[i];
        path_length += graph.get_length(handle);
        handles_iv[i] = as_integer(local_handle(handle));
        // we will explode if the node isn't in the graph
    }
    sdsl::util::assign(handles, sdsl::enc_vector<>(handles_iv));

    // make the bitvector for path offsets
    sdsl::bit_vector offsets_bv;
    sdsl::util::assign(offsets_bv, sdsl::bit_vector(path_length));

    //cerr << "path " << path_name << " has " << path.size() << endl;
    for (size_t i = 0; i < path.size(); ++i) {
        //cerr << i << endl;
        auto& handle = path[i];
        // record position of node
        offsets_bv[path_off] = 1;
        // and update the offset counter
        path_off += graph.get_length(handle);
    }
    // and path offsets
    sdsl::util::assign(offsets, sdsl::rrr_vector<>(offsets_bv));
    // and set up rank/select dictionary on them
    sdsl::util::assign(offsets_rank, sdsl::rrr_vector<>::rank_1_type(&offsets));
    sdsl::util::assign(offsets_select, sdsl::rrr_vector<>::select_1_type(&offsets));
}

handle_t XGPath::handle(size_t offset) const {
    return external_handle(as_handle(handles[offset]));
}
    
handle_t XGPath::handle_at_position(size_t pos) const {
    return handle(step_rank_at_position(pos));
}

size_t XGPath::handle_start(size_t offset) const {
    return offsets_select(offset+1);
}

size_t XGPath::step_rank_at_position(size_t pos) const {
    return offsets_rank(pos+1)-1;
}

bool XGPath::is_reverse(size_t offset) const {
    return number_bool_packing::unpack_bit(as_handle(handles[offset]));
}

handle_t XGPath::local_handle(const handle_t& handle) const {
    if (as_integer(handle) < as_integer(min_handle)) {
        return as_handle(std::numeric_limits<int64_t>::max());
    } else {
        return as_handle(as_integer(handle)-as_integer(min_handle));
    }
}

handle_t XGPath::external_handle(const handle_t& handle) const {
    return as_handle(as_integer(handle)+as_integer(min_handle));
}

void XG::serialize(ostream& out) const {
    serialize_and_measure(out);
}

size_t XG::serialize_and_measure(ostream& out, sdsl::structure_tree_node* s, std::string name) const {

    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
    size_t written = 0;
    
    // Do the magic number
    out << "XG";
    written += 2;
    
    // And the version number
    int32_t version_buffer = htonl(CURRENT_VERSION);
    out.write((char*) &version_buffer, sizeof(version_buffer));
    written += sizeof(version_buffer) / sizeof(char);
    
    ////////////////////////////////////////////////////////////////////////
    // DO NOT CHANGE THIS CODE without creating a new XG version:
    // 1. Increment CURRENT_VERSION to a new integer.
    // 2. Add your new serialization code.
    // 3. Add a new case for your new version to XG::load()
    ////////////////////////////////////////////////////////////////////////

    written += sdsl::write_member(s_iv.size(), out, child, "sequence_length");
    written += sdsl::write_member(node_count, out, child, "node_count");
    written += sdsl::write_member(edge_count, out, child, "edge_count");
    written += sdsl::write_member(path_count, out, child, "path_count");
    written += sdsl::write_member(min_id, out, child, "min_id");
    written += sdsl::write_member(max_id, out, child, "max_id");

    written += r_iv.serialize(out, child, "rank_id_vector");

    written += g_iv.serialize(out, child, "graph_vector");
    written += g_bv.serialize(out, child, "graph_bit_vector");
    written += g_bv_rank.serialize(out, child, "graph_bit_vector_rank");
    written += g_bv_select.serialize(out, child, "graph_bit_vector_select");
    
    written += s_iv.serialize(out, child, "seq_vector");
    written += s_bv.serialize(out, child, "seq_node_starts");
    written += s_bv_rank.serialize(out, child, "seq_node_starts_rank");
    written += s_bv_select.serialize(out, child, "seq_node_starts_select");

    // Treat the paths as their own node
    size_t paths_written = 0;
    auto paths_child = sdsl::structure_tree::add_child(child, "paths", sdsl::util::class_name(*this));

    paths_written += pn_iv.serialize(out, paths_child, "path_names");
    paths_written += pn_csa.serialize(out, paths_child, "path_names_csa");
    paths_written += pn_bv.serialize(out, paths_child, "path_names_starts");
    paths_written += pn_bv_rank.serialize(out, paths_child, "path_names_starts_rank");
    paths_written += pn_bv_select.serialize(out, paths_child, "path_names_starts_select");
    paths_written += pi_iv.serialize(out, paths_child, "path_ids");
    // TODO: Path count is written twice (once from paths.size() and once earlier from path_count)
    // We should remove one and cut a new xg version
    paths_written += sdsl::write_member(paths.size(), out, paths_child, "path_count");    
    for (size_t i = 0; i < paths.size(); i++) {
        XGPath* path = paths[i];
        paths_written += path->serialize(out, paths_child, "path:" + get_path_name(as_path_handle(i+1)));
    }

    paths_written += np_bv.serialize(out, paths_child, "node_path_mapping_starts");
    paths_written += np_bv_select.serialize(out, paths_child, "node_path_mapping_starts_select");
    paths_written += np_iv.serialize(out, paths_child, "node_path_mapping");
    paths_written += nr_iv.serialize(out, paths_child, "node_path_rank");
    paths_written += nx_iv.serialize(out, paths_child, "node_path_position");
    
    sdsl::structure_tree::add_size(paths_child, paths_written);
    written += paths_written;

    sdsl::structure_tree::add_size(child, written);
    return written;
    
}

void XG::dump_to_stream(std::ostream& out) const {
    out << "===BEGIN XG===" << std::endl;
    out << "index\tg\tbv" << std::endl;
    for (size_t i = 0; i < g_iv.size(); i++) {
        out << i << "\t" << g_iv[i] << "\t" << g_bv[i] << std::endl;
    }
    out << "===END XG===" << std::endl;
}

void XG::from_gfa(const std::string& gfa_filename, bool validate, std::string basename) {
    if (basename.empty()) {
        basename = temp_file::create();
    }
    char* filename = (char*)gfa_filename.c_str();
    gfak::GFAKluge gfa;
    // set up our enumerators
    auto for_each_sequence = [&](const std::function<void(const std::string& seq, const nid_t& node_id)>& lambda) {
        gfa.for_each_sequence_line_in_file(filename, [&](gfak::sequence_elem s) {
            nid_t node_id = std::stol(s.name);
            lambda(s.sequence, node_id);
        });
    };
    auto for_each_edge = [&](const std::function<void(const nid_t& from_id, const bool& from_rev,
                                                      const nid_t& to_id, const bool& to_rev)>& lambda) {
        gfa.for_each_edge_line_in_file(filename, [&](gfak::edge_elem e) {
            if (e.source_name.empty()) return;
            nid_t from_id = std::stol(e.source_name);
            bool from_rev = !e.source_orientation_forward;
            nid_t to_id = std::stol(e.sink_name);
            bool to_rev = !e.sink_orientation_forward;
            lambda(from_id, from_rev, to_id, to_rev);
        });
    };
    auto for_each_path_element = [&](const std::function<void(const std::string& path_name,
                                                              const nid_t& node_id, const bool& is_rev,
                                                              const std::string& cigar,
                                                              const bool& is_empty, const bool& is_circular)>& lambda) {
        gfa.for_each_path_element_in_file(filename, [&](const std::string& path_name_raw, const std::string& node_id_str,
                                                        bool is_rev, const std::string& cigar,
                                                        bool is_empty, bool is_circular) {
            nid_t node_id = std::stol(node_id_str);
            std::string path_name = path_name_raw;
            path_name.erase(std::remove_if(path_name.begin(), path_name.end(), [](char c) { return std::isspace(c); }), path_name.end());
            lambda(path_name, node_id, is_rev, cigar, is_empty, is_circular);
        });
    };
    from_enumerators(for_each_sequence, for_each_edge, for_each_path_element, validate);
}


/// build the graph from another simple graph
void XG::from_handle_graph(const HandleGraph& graph) {
    // set up our enumerators
    auto for_each_sequence = [&](const std::function<void(const std::string& seq, const nid_t& node_id)>& lambda) {
        graph.for_each_handle([&](const handle_t& handle) {
                lambda(graph.get_sequence(handle), graph.get_id(handle));
            });
    };
    auto for_each_edge = [&](const std::function<void(const nid_t& from_id, const bool& from_rev,
                                                      const nid_t& to_id, const bool& to_rev)>& lambda) {
        graph.for_each_edge([&](const edge_t& edge) {
                lambda(graph.get_id(edge.first), graph.get_is_reverse(edge.first),
                       graph.get_id(edge.second), graph.get_is_reverse(edge.second));
            });
    };
    auto for_each_path_element = [&](const std::function<void(const std::string& path_name,
                                                              const nid_t& node_id, const bool& is_rev,
                                                              const std::string& cigar, const bool& is_empty, const bool& is_circular)>& lambda) {
        // no-op
    };
    from_enumerators(for_each_sequence, for_each_edge, for_each_path_element, false);
}

/// build the graph from another path handle graph
void XG::from_path_handle_graph(const PathHandleGraph& graph) {
    // set up our enumerators
    auto for_each_sequence = [&](const std::function<void(const std::string& seq, const nid_t& node_id)>& lambda) {
        graph.for_each_handle([&](const handle_t& handle) {
                lambda(graph.get_sequence(handle), graph.get_id(handle));
            });
    };
    auto for_each_edge = [&](const std::function<void(const nid_t& from_id, const bool& from_rev,
                                                      const nid_t& to_id, const bool& to_rev)>& lambda) {
        graph.for_each_edge([&](const edge_t& edge) {
                lambda(graph.get_id(edge.first), graph.get_is_reverse(edge.first),
                       graph.get_id(edge.second), graph.get_is_reverse(edge.second));
            });
    };
    auto for_each_path_element = [&](const std::function<void(const std::string& path_name,
                                                              const nid_t& node_id, const bool& is_rev,
                                                              const std::string& cigar, const bool& is_empty, const bool& is_circular)>& lambda) {
        graph.for_each_path_handle([&](const path_handle_t& path_handle) {
                std::string path_name = graph.get_path_name(path_handle);
                size_t step_count = 0;
                bool path_is_circular = graph.get_is_circular(path_handle);
                graph.for_each_step_in_path(path_handle, [&](const step_handle_t& step) {
                        handle_t handle = graph.get_handle_of_step(step);
                        lambda(path_name, graph.get_id(handle), graph.get_is_reverse(handle), "", false, path_is_circular);
                        ++step_count;
                    });
                if (step_count == 0) {
                    lambda(path_name, 0, false, "", true, path_is_circular);
                }
            });
    };
    from_enumerators(for_each_sequence, for_each_edge, for_each_path_element, false);
}

void XG::from_enumerators(const std::function<void(const std::function<void(const std::string& seq, const nid_t& node_id)>&)>& for_each_sequence,
                          const std::function<void(const std::function<void(const nid_t& from, const bool& from_rev,
                                                                            const nid_t& to, const bool& to_rev)>&)>& for_each_edge,
                          const std::function<void(const std::function<void(const std::string& path_name,
                                                                            const nid_t& node_id, const bool& is_rev,
                                                                            const std::string& cigar, const bool& is_empty,
                                                                            const bool& is_circular)>&)>& for_each_path_element,
                          bool validate, std::string basename) {

    if (basename.empty()) {
        basename = temp_file::create();
    }
    node_count = 0;
    seq_length = 0;
    edge_count = 0;
    path_count = 0;
    min_id = std::numeric_limits<int64_t>::max();
    max_id = 0;
    // get information about graph size and id ranges
#ifdef VERBOSE_DEBUG
    cerr << "computing graph sequence length and node count" << endl;
#endif
    for_each_sequence([&](const std::string& seq, const nid_t& id) {
            // min id starts at 0
            min_id = std::min(min_id, id);
            max_id = std::max(max_id, id);
            seq_length += seq.size();
            ++node_count;
        });
#ifdef VERBOSE_DEBUG
    cerr << "counting edges" << endl;
#endif
    // edge count
    for_each_edge([&](const nid_t& from_id, const bool& from_rev, const nid_t& to_id, const bool& to_rev) {
            ++edge_count;
        });
    // path count
    std::string pname;
#ifdef VERBOSE_DEBUG
    cerr << "counting paths" << endl;
#endif
    for_each_path_element([&](const std::string& path_name, const nid_t& node_id, const bool& is_rev, const std::string& cigar, const bool& is_empty, const bool& is_circular) {
            if (path_name != pname) {
                ++path_count;
            }
            pname = path_name;
        });
#ifdef VERBOSE_DEBUG
    std::cerr << "graph has " << seq_length << "bp in sequence, "
              << node_count << " nodes, "
              << edge_count << " edges, and "
              << path_count << " paths " << std::endl;
#endif

    // set up our compressed representation
    sdsl::int_vector<> i_iv;
    sdsl::util::assign(s_iv, sdsl::int_vector<>(seq_length, 0, 3));
    sdsl::util::assign(s_bv, sdsl::bit_vector(seq_length+1));
    sdsl::util::assign(i_iv, sdsl::int_vector<>(node_count));
    sdsl::util::assign(r_iv, sdsl::int_vector<>(max_id-min_id+1)); // note: possibly discontiguous
    
    // for each node in the sequence
    // concatenate the labels into the s_iv
#ifdef VERBOSE_DEBUG
    cerr << "storing node labels" << endl;
#endif
    size_t r = 1;    
    // first make i_iv and r_iv
    for_each_sequence([&](const std::string& seq, const nid_t& id) {
            i_iv[r-1] = id;
            // store ids to rank mapping
            r_iv[id-min_id] = r;
            ++r;
        });
    sdsl::util::bit_compress(i_iv);
    sdsl::util::bit_compress(r_iv);

    // then make s_bv and s_iv
    size_t j = 0;
    for_each_sequence([&](const std::string& seq, const nid_t& id) {
            //size_t i = r_iv[id-min_id]-1;
            s_bv[j] = 1; // record node start
            for (auto c : seq) {
                s_iv[j++] = dna3bit(c); // store sequence
            }
        });
    s_bv[seq_length] = 1;

    // to label the paths we'll need to compress and index our vectors
    sdsl::util::bit_compress(s_iv);
    sdsl::util::assign(s_bv_rank, sdsl::rank_support_v<1>(&s_bv));
    sdsl::util::assign(s_bv_select, sdsl::bit_vector::select_1_type(&s_bv));
    
    // now that we've set up our sequence indexes, we can build the locally traversable graph storage

    auto temp_get_handle = [&](const nid_t& id, bool orientation) {
        uint64_t handle_rank = r_iv[id-min_id];
        return number_bool_packing::pack(handle_rank, orientation);
    };
    auto temp_node_size = [&](const nid_t& id) {
        uint64_t handle_rank = r_iv[id-min_id];
        return s_bv_select(handle_rank+1)-s_bv_select(handle_rank);
    };
    auto temp_get_id = [&](const handle_t& h) {
        return i_iv[number_bool_packing::unpack_number(h)-1];
    };

#ifdef VERBOSE_DEBUG
    cerr << "collecting edges " << endl;
#endif
    
    // first, we need to collect the edges for each node
    // we use the mmmultimap here to reduce in-memory costs to a minimum
    std::string edge_f_t_idx = basename + ".from_to.mm";
    std::string edge_t_f_idx = basename + ".to_from.mm";
    auto edge_from_to_mm = std::make_unique<mmmulti::map<uint64_t, uint64_t>>(edge_f_t_idx);
    auto edge_to_from_mm = std::make_unique<mmmulti::map<uint64_t, uint64_t>>(edge_t_f_idx);
    for_each_edge([&](const nid_t& from_id, const bool& from_rev, const nid_t& to_id, const bool& to_rev) {
            handle_t from_handle = temp_get_handle(from_id, from_rev);
            handle_t to_handle = temp_get_handle(to_id, to_rev);
            edge_from_to_mm->append(as_integer(from_handle), as_integer(to_handle));
            edge_to_from_mm->append(as_integer(to_handle), as_integer(from_handle));
        });
    handle_t max_handle = number_bool_packing::pack(r_iv.size(), true);
    edge_from_to_mm->index(as_integer(max_handle));
    edge_to_from_mm->index(as_integer(max_handle));

    // calculate g_iv size
    size_t g_iv_size =
        node_count * G_NODE_HEADER_LENGTH // record headers
        + edge_count * 2 * G_EDGE_LENGTH; // edges (stored twice)
    sdsl::util::assign(g_iv, sdsl::int_vector<>(g_iv_size));
    sdsl::util::assign(g_bv, sdsl::bit_vector(g_iv_size));

#ifdef VERBOSE_DEBUG
    cerr << "computing graph vector " << endl;
#endif
    
    int64_t g = 0; // pointer into g_iv and g_bv
    for (int64_t i = 0; i < node_count; ++i) {
        nid_t id = i_iv[i];
#ifdef VERBOSE_DEBUG
        if (i % 1000 == 0) cerr << i << " of " << node_count << " ~ " << (float)i/(float)node_count * 100 << "%" << "\r";
#endif
        handle_t handle = temp_get_handle(id, false);
        //std::cerr << "id " << id << std::endl;
        g_bv[g] = 1; // mark record start for later query
        g_iv[g++] = id;
        g_iv[g++] = node_vector_offset(id);
        g_iv[g++] = temp_node_size(id);
        size_t to_edge_count = 0;
        size_t from_edge_count = 0;
        size_t to_edge_count_idx = g++;
        size_t from_edge_count_idx = g++;
        // write the edges in id-based format
        // we will next convert these into relative format
        for (auto orientation : { false, true }) {
            handle_t to = temp_get_handle(id, orientation);
            //std::cerr << "looking at to handle " << number_bool_packing::unpack_number(to) << ":" << number_bool_packing::unpack_bit(to) << std::endl;
            edge_to_from_mm->for_unique_values_of(as_integer(to), [&](const uint64_t& _from) {
                    handle_t from = as_handle(_from);
                    //std::cerr << "edge to " << number_bool_packing::unpack_number(from) << ":" << number_bool_packing::unpack_bit(from)
                    //<< " -> " << number_bool_packing::unpack_number(to) << ":" << number_bool_packing::unpack_bit(to) << std::endl;
                    g_iv[g++] = temp_get_id(from);
                    g_iv[g++] = edge_type(from, to);
                    ++to_edge_count;
                });
        }
        g_iv[to_edge_count_idx] = to_edge_count;
        for (auto orientation : { false, true }) {
            handle_t from = temp_get_handle(id, orientation);
            //std::cerr << "looking at from handle " << number_bool_packing::unpack_number(from) << ":" << number_bool_packing::unpack_bit(from) << std::endl;
            edge_from_to_mm->for_unique_values_of(as_integer(from), [&](const uint64_t& _to) {
                    handle_t to = as_handle(_to);
                    //std::cerr << "edge from " << number_bool_packing::unpack_number(from) << ":" << number_bool_packing::unpack_bit(from)
                    //<< " -> " << number_bool_packing::unpack_number(to) << ":" << number_bool_packing::unpack_bit(to) << std::endl;
                    g_iv[g++] = temp_get_id(to);
                    g_iv[g++] = edge_type(from, to);
                    ++from_edge_count;
                });
        }
        g_iv[from_edge_count_idx] = from_edge_count;
    }

#ifdef VERBOSE_DEBUG
    std::cerr << node_count << " of " << node_count << " ~ 100.0000%" << std::endl;
#endif

    // cleanup our mmmultimap
    edge_from_to_mm.reset();
    edge_to_from_mm.reset();
    std::remove(edge_f_t_idx.c_str());
    std::remove(edge_t_f_idx.c_str());

    // set up rank and select supports on g_bv so we can locate nodes in g_iv
    sdsl::util::assign(g_bv_rank, sdsl::rank_support_v<1>(&g_bv));
    sdsl::util::assign(g_bv_select, sdsl::bit_vector::select_1_type(&g_bv));

#ifdef VERBOSE_DEBUG
    cerr << "making graph vector relativistic " << endl;
#endif

    // convert the edges in g_iv to relativistic form
    for (int64_t i = 0; i < node_count; ++i) {
        int64_t id = i_iv[i];
        // find the start of the node's record in g_iv
        int64_t g = g_bv_select(id_to_rank(id));
        // get to the edges to
        int edges_to_count = g_iv[g+G_NODE_TO_COUNT_OFFSET];
        int edges_from_count = g_iv[g+G_NODE_FROM_COUNT_OFFSET];
        int64_t t = g + G_NODE_HEADER_LENGTH;
        int64_t f = g + G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * edges_to_count;
        for (int64_t j = t; j < f; ) {
            g_iv[j] = g_bv_select(id_to_rank(g_iv[j])) - g;
            j += 2;
        }
        for (int64_t j = f; j < f + G_EDGE_LENGTH * edges_from_count; ) {
            g_iv[j] = g_bv_select(id_to_rank(g_iv[j])) - g;
            j += 2;
        }
    }
    sdsl::util::clear(i_iv);
    sdsl::util::bit_compress(g_iv);

#ifdef VERBOSE_DEBUG
    cerr << "storing paths" << endl;
#endif
    // paths
    std::string path_names;

    std::string curr_path_name;
    std::vector<handle_t> curr_path_steps;
    size_t curr_node_count = 0;
    bool curr_is_circular = false; // TODO, use TP:Z:circular tag... we'll have to fish this out of the file
    uint64_t p = 0;

    auto build_accumulated_path = [&](void) {
        // only build if we had a path to build
#ifdef VERBOSE_DEBUG
        if (++p % 100 == 0) std::cerr << p << " of " << path_count << " ~ " << (float)p/(float)path_count * 100 << "%" << "\r";
#endif
        size_t unique_member_count = 0;
        path_names += start_marker + curr_path_name + end_marker;
        XGPath* path = new XGPath(curr_path_name, curr_path_steps,
                                  curr_is_circular,
                                  *this);
        paths.push_back(path);
    };

    // todo ... is it circular?
    // might make sense to scan the file for this
    bool has_path = false;
    for_each_path_element([&](const std::string& path_name, const nid_t& node_id, const bool& is_rev, const std::string& cigar, const bool& is_empty, const bool& is_circular) {
            if (path_name != curr_path_name && !curr_path_name.empty()) {
                // build the last path we've accumulated
                build_accumulated_path();
                curr_path_steps.clear();
                curr_is_circular = false;
            }
            curr_path_name = path_name;
            if (!is_empty) {
                handle_t visiting = get_handle(node_id, is_rev);
#ifdef debug_path_index
                std::cerr << "Adding handle for " << node_id << (is_rev ? "-" : "+") << " with value "
                    << as_integer(visiting) << " to path " << path_name << std::endl;
                std::cerr << "Node length is " << get_length(visiting) << " bp" << std::endl;
#endif
                curr_path_steps.push_back(visiting);
            }
            curr_is_circular = is_circular;
            has_path = true;
        });
    // build the last path
    if (has_path) {
        build_accumulated_path();
    }
    curr_path_steps.clear();
    curr_is_circular = false;
#ifdef VERBOSE_DEBUG
    std::cerr << path_count << " of " << path_count << " ~ 100.0000%" << std::endl;
#endif

    //std::cerr << "path names " << path_names << std::endl;

    // handle path names
    sdsl::util::assign(pn_iv, sdsl::int_vector<>(path_names.size()));
    sdsl::util::assign(pn_bv, sdsl::bit_vector(path_names.size()));
    // now record path name starts
    for (size_t i = 0; i < path_names.size(); ++i) {
        pn_iv[i] = path_names[i];
        if (path_names[i] == start_marker) {
            pn_bv[i] = 1; // register name start
        }
    }
    sdsl::util::assign(pn_bv_rank, sdsl::rank_support_v<1>(&pn_bv));
    sdsl::util::assign(pn_bv_select, sdsl::bit_vector::select_1_type(&pn_bv));
    
    // is this file removed by construct?
    string path_name_file = basename + ".pathnames.iv";
    sdsl::store_to_file((const char*)path_names.c_str(), path_name_file);
    sdsl::construct(pn_csa, path_name_file, 1);

#ifdef VERBOSE_DEBUG
    cerr << "computing node to path membership" << endl;
#endif
    
    // create the node-to-path indexes
    index_node_to_path(basename);

    // validate the graph
    if (validate) {
        // do we get the correct sequences when looking up handles by id?
        for_each_sequence([&](const std::string& seq, const nid_t& id) {
                handle_t handle = get_handle(id);
                if (seq != get_sequence(handle)) {
                    std::cerr << "mismatch in handle sequence for " << id << std::endl;
                    exit(1);
                }
                if (get_id(handle) != id) {
                    std::cerr << "mismatch in id for " << id << std::endl;
                    exit(1);
                }
            });
        // do we have the correct set of edges?
        for_each_edge([&](const nid_t& from_id, const bool& from_rev, const nid_t& to_id, const bool& to_rev) {
                handle_t from_handle = get_handle(from_id, from_rev);
                handle_t to_handle = get_handle(to_id, to_rev);
                bool seen_to = false;
                follow_edges(from_handle, false, [&](const handle_t& h) {
                        //std::cerr << "fwd I see edge " << get_id(from_handle) << ":" << get_is_reverse(from_handle) << " -> " << get_id(h) << ":" << get_is_reverse(h) << std::endl;
                        //std::cerr << as_integer(h) << " ==? " << as_integer(to_handle) << std::endl;
                        if (h == to_handle) {
                            seen_to = true;
                        }
                    });
                bool seen_from = false;
                follow_edges(to_handle, true, [&](const handle_t& h) {
                        //std::cerr << "rev I see edge " << get_id(h) << ":" << get_is_reverse(h) << " -> " << get_id(to_handle) << ":" << get_is_reverse(to_handle) << std::endl;
                        //std::cerr << as_integer(h) << " ==? " << as_integer(from_handle) << std::endl;
                        if (h == from_handle) {
                            seen_from = true;
                        }
                    });
                if (!seen_to) {
                    std::cerr << "can't find to edge for " << get_id(from_handle) << ":" << get_is_reverse(from_handle)
                              << " -> " << get_id(to_handle) << ":" << get_is_reverse(to_handle) << std::endl;
                    exit(1);
                }
                if (!seen_from) {
                    std::cerr << "can't find from edge for " << get_id(from_handle) << ":" << get_is_reverse(from_handle)
                              << " -> " << get_id(to_handle) << ":" << get_is_reverse(to_handle) << std::endl;
                    exit(1);
                }
            });
        // do our stored paths match those in the input?

        std::string curr_path_name;
        std::vector<handle_t> curr_path_steps;
        size_t curr_node_count = 0;
        bool curr_is_circular = false; // TODO, use TP:Z:circular tag... we'll have to fish this out of the file
        uint64_t p_handle = 0;
        auto check_accumulated_path = [&](void) {
            ++p_handle;
            // only check if we had a path to build
            if (curr_path_steps.empty()) return;
            auto& path = *paths[p_handle-1];
            // check that the path name is correct
            if (get_path_name(as_path_handle(p_handle)) != curr_path_name) {
                std::cerr << "path name mismatch " << get_path_name(as_path_handle(p_handle)) << " != " << curr_path_name << std::endl;
                exit(1);
            }
            // check that the path handles are correct
            for (uint64_t i = 0; i < curr_path_steps.size(); ++i) {
                if (path.handle(i) != curr_path_steps[i]) {
                    std::cerr << "handle mismatch " << get_id(path.handle(i)) << " != " << get_id(curr_path_steps[i]) << " in path " << curr_path_name << std::endl;
                    exit(1);
                }
            }
            // check that the path positions are correct
            uint64_t pos = 0;
            path_handle_t path_handle = as_path_handle(p_handle);
            //std::cerr << "on path " << curr_path_name << std::endl;
            //const std::function<bool(const step_handle_t&, const bool&, const uint64_t&)>& iteratee) const;
            for (uint64_t i = 0; i < curr_path_steps.size(); ++i) {
                handle_t handle = curr_path_steps[i];
                //std::cerr << "looking at node " << get_id(handle) << " on path " << curr_path_name << std::endl;
                uint64_t handle_length = get_length(handle);
                for (uint64_t j = 0; j < handle_length; ++j) {
                    handle_t handle_at = path.handle_at_position(pos+j);
                    if (handle != handle_at) {
                        std::cerr << "handle at position mismatch " << get_id(handle) << " != " << get_id(handle_at)
                                  << " in path " << curr_path_name << " at position " << pos+j << std::endl;
                        exit(1);
                    }
                }
                bool path_seen = false;
                auto check_pos_index = [&](const step_handle_t& step, const bool& is_rev, const uint64_t& pos) {
                    // check that the step handle is the same as this handle
                    //std::cerr << "checking " << get_id(handle) << " on " << curr_path_name << std::endl;
                    path_handle_t path_handle_of = get_path_handle_of_step(step);
                    if (path_handle == path_handle_of && as_integers(step)[1] == i) {
                        //std::cerr << "found matching step handle" << std::endl;
                        handle_t oriented_handle = !get_is_reverse(handle) && is_rev ? flip(handle) : handle;
                        handle_t handle_at = get_handle_of_step(step);
                        if (oriented_handle != handle_at) {
                            std::cerr << "handle at step mismatch " << get_id(oriented_handle) << " != " << get_id(handle_at)
                                      << " in path " << curr_path_name << std::endl;
                            std::cerr << "handle " << as_integer(handle_at) << " != " << as_integer(oriented_handle) << std::endl;
                            exit(1);
                        }
                        handle_t handle_at_pos = path.handle_at_position(pos);
                        if (handle_at != handle_at_pos) {
                            std::cerr << "handle at position mismatch " << get_id(handle_at) << " != " << get_id(handle_at_pos)
                                      << " in path " << curr_path_name << " at position " << pos << std::endl;
                            exit(1);
                        }
                        path_seen = true;
                        return false;
                    }
                    return true;
                };
                // verify that we have the right position in the node to path position index
                for_each_step_position_on_handle(handle, check_pos_index);
                if (!path_seen) {
                    std::cerr << "didn't find path " << curr_path_name << " in reverse index for " << get_id(handle) << std::endl;
                    exit(1);
                }
                pos += get_length(handle);
            }
        };
        for_each_path_element([&](const std::string& path_name, const nid_t& node_id, const bool& is_rev, const std::string& cigar, const bool& is_empty, const bool& is_circular) {
                if (path_name != curr_path_name && !curr_path_name.empty()) {
                    // check the last path we've accumulated
                    check_accumulated_path();
                    curr_path_steps.clear();
                }
                curr_path_name = path_name;
                if (!is_empty) {
                    curr_path_steps.push_back(get_handle(node_id, is_rev));
                }
            });
        // check the last path
        check_accumulated_path();
        curr_path_steps.clear();
    }

//#define DEBUG_CONSTRUCTION
#ifdef DEBUG_CONSTRUCTION
    cerr << "|g_iv| = " << size_in_mega_bytes(g_iv) << endl;
    cerr << "|g_bv| = " << size_in_mega_bytes(g_bv) << endl;
    cerr << "|s_iv| = " << size_in_mega_bytes(s_iv) << endl;

    cerr << "|r_iv| = " << size_in_mega_bytes(r_iv) << endl;

    cerr << "|s_bv| = " << size_in_mega_bytes(s_bv) << endl;
    
    long double paths_mb_size = 0;
    cerr << "|pn_iv| = " << size_in_mega_bytes(pn_iv) << endl;
    paths_mb_size += size_in_mega_bytes(pn_iv);
    cerr << "|pn_csa| = " << size_in_mega_bytes(pn_csa) << endl;
    paths_mb_size += size_in_mega_bytes(pn_csa);
    cerr << "|pn_bv| = " << size_in_mega_bytes(pn_bv) << endl;
    paths_mb_size += size_in_mega_bytes(pn_bv);
    paths_mb_size += size_in_mega_bytes(pn_bv_rank);
    paths_mb_size += size_in_mega_bytes(pn_bv_select);
    paths_mb_size += size_in_mega_bytes(pi_iv);
    cerr << "|np_iv| = " << size_in_mega_bytes(np_iv) << endl;
    paths_mb_size += size_in_mega_bytes(np_iv);
    cerr << "|np_bv| = " << size_in_mega_bytes(np_bv) << endl;
    paths_mb_size += size_in_mega_bytes(np_bv);
    paths_mb_size += size_in_mega_bytes(np_bv_select);
    cerr << "total paths size " << paths_mb_size << endl;

    float path_ids_mb_size=0;
    float path_dir_mb_size=0;
    float path_pos_mb_size=0;
    float path_ranks_mb_size=0;
    float path_offsets_mb_size=0;
    for (size_t i = 0; i < paths.size(); i++) {
        // Go through paths by number, so we can determine rank
        XGPath* path = paths[i];
        path_ids_mb_size += size_in_mega_bytes(path->ids);
        path_dir_mb_size += size_in_mega_bytes(path->directions);
        path_pos_mb_size += size_in_mega_bytes(path->positions);
        path_ranks_mb_size += size_in_mega_bytes(path->ranks);
        path_offsets_mb_size += size_in_mega_bytes(path->offsets);
    }
    cerr << "path ids size " << path_ids_mb_size << endl;
    cerr << "path directions size " << path_dir_mb_size << endl;
    cerr << "path positions size " << path_pos_mb_size << endl;
    cerr << "path ranks size " << path_ranks_mb_size << endl;
    cerr << "path offsets size " << path_offsets_mb_size << endl;

#endif

    bool print_graph = false;
    if (print_graph) {
        cerr << "printing graph" << endl;
        // we have to print the relativistic graph manually because the default sdsl printer assumes unsigned integers are stored in it
        for (size_t i = 0; i < g_iv.size(); ++i) {
            cerr << (int64_t)g_iv[i] << " ";
        } cerr << endl;
        for (int64_t i = 0; i < node_count; ++i) {
            int64_t id = rank_to_id(i+1);
            // find the start of the node's record in g_iv
            int64_t g = g_bv_select(id_to_rank(id));
            // get to the edges to
            int edges_to_count = g_iv[g+G_NODE_TO_COUNT_OFFSET];
            int edges_from_count = g_iv[g+G_NODE_FROM_COUNT_OFFSET];
            int sequence_size = g_iv[g+G_NODE_LENGTH_OFFSET];
            size_t seq_start = g_iv[g+G_NODE_SEQ_START_OFFSET];
            cerr << id << " ";
            for (int64_t j = seq_start; j < seq_start+sequence_size; ++j) {
                cerr << revdna3bit(s_iv[j]);
            } cerr << " : ";
            int64_t t = g + G_NODE_HEADER_LENGTH;
            int64_t f = g + G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * edges_to_count;
            cerr << " from ";
            for (int64_t j = t; j < f; ) {
                cerr << rank_to_id(g_bv_rank(g+g_iv[j])+1) << " ";
                j += 2;
            }
            cerr << " to ";
            for (int64_t j = f; j < f + G_EDGE_LENGTH * edges_from_count; ) {
                cerr << rank_to_id(g_bv_rank(g+g_iv[j])+1) << " ";
                j += 2;
            }
            cerr << endl;
        }
        cerr << s_iv << endl;
        for (size_t i = 0; i < s_iv.size(); ++i) {
            cerr << revdna3bit(s_iv[i]);
        } cerr << endl;
        cerr << s_bv << endl;
        cerr << "paths (" << paths.size() << ")" << endl;
        for (size_t i = 0; i < paths.size(); i++) {
            // Go through paths by number, so we can determine rank
            XGPath* path = paths[i];
            cerr << get_path_name(as_path_handle(i + 1)) << endl;
            // manually print IDs because simplified wavelet tree doesn't support ostream for some reason
            for (size_t j = 0; j + 1 < path->handles.size(); j++) {
                cerr << get_id(path->handle(j)) << " ";
            }
            if (path->handles.size() > 0) {
                cerr << get_id(path->handle(path->handles.size() - 1));
            }
            cerr << endl;
            cerr << path->offsets << endl;
        }
        cerr << np_bv << endl;
        cerr << np_iv << endl;
        cerr << nx_iv << endl;
        
    }

}

void XG::index_node_to_path(const std::string& basename) {
    
    // node -> paths
    // use the mmmultimap...
    std::string node_path_idx = basename + ".node_path.mm";
    auto node_path_mm = std::make_unique<mmmulti::map<uint64_t, std::tuple<uint64_t, uint64_t, uint64_t>>>(node_path_idx);
    uint64_t path_step_count = 0;
    // for each path...
    // could be done in parallel
    for (size_t i = 1; i <= paths.size(); ++i) {
        // Go through paths by number, so we can determine rank
        path_handle_t path_handle = as_path_handle(i);
        const XGPath& path = *paths[i-1];
#ifdef debug_path_index
        std::cerr << "Indexing path " << &path << " at index " << i-1 << " of " << paths.size() << std::endl; 
#endif
        uint64_t pos = 0;
#ifdef VERBOSE_DEBUG
        if (i % 100 == 0) std::cerr << i << " of " << path_count << " ~ " << (float)i/(float)path_count * 100 << "%" << "\r";
#endif
        for (size_t j = 0; j < path.handles.size(); ++j) {
            handle_t handle = path.handle(j);
            
#ifdef debug_path_index
            std::cerr << "Look at handle " << j << " of " << path.handles.size() << ": "
                << get_id(handle) << (get_is_reverse(handle) ? "-" : "+") << " with value " << as_integer(handle) << std::endl;
                
            std::cerr << "Handle is for node at g vector position " << handlegraph::number_bool_packing::unpack_number(handle)
                << " of " << g_iv.size() << " and node length will be at "
                << handlegraph::number_bool_packing::unpack_number(handle) + G_NODE_LENGTH_OFFSET << std::endl;
#endif
            uint64_t handle_length = get_length(handle);
            bool is_rev = number_bool_packing::unpack_bit(handle);
            // and record the relative orientation by packing it into the position
            uint64_t path_and_rev = as_integer(number_bool_packing::pack(as_integer(path_handle), is_rev));
            // determine the path relative position on the forward strand
            uint64_t adj_pos = is_rev ? pos + handle_length - 1: pos;
            node_path_mm->append(id_to_rank(get_id(handle)),
                                std::make_tuple(path_and_rev, j, adj_pos));
            ++path_step_count;
            pos += handle_length;
        }
    }
#ifdef VERBOSE_DEBUG
    std::cerr << path_count << " of " << path_count << " ~ 100.0000%" << std::endl;
#endif
    node_path_mm->index(node_count+1);
    
#ifdef VERBOSE_DEBUG
    std::cerr << "determining size of node to path position mappings" << std::endl;
#endif
    uint64_t np_size = 0;
    for (int64_t i = 0; i < node_count; ++i) {
#ifdef VERBOSE_DEBUG
        if (i % 1000 == 0) cerr << i << " of " << node_count << " ~ " << (float)i/(float)node_count * 100 << "%" << "\r";
#endif
        uint64_t has_steps = false;
        node_path_mm->for_values_of(i+1, [&](const std::tuple<uint64_t, uint64_t, uint64_t>& v) {
            ++np_size;
            has_steps = true;
        });
        if (!has_steps) ++np_size; // skip over null entry
    }
#ifdef VERBOSE_DEBUG
    std::cerr << node_count << " of " << node_count << " ~ 100.0000%" << std::endl;
#endif
    sdsl::util::assign(np_iv, sdsl::int_vector<>(np_size));
    sdsl::util::assign(nr_iv, sdsl::int_vector<>(np_size));
    sdsl::util::assign(nx_iv, sdsl::int_vector<>(np_size));
    sdsl::util::assign(np_bv, sdsl::bit_vector(np_size));
    // now we need to map this into a new data structure for node->path position mappings
#ifdef VERBOSE_DEBUG
    std::cerr << "recording node to path position mappings" << std::endl;
#endif
    uint64_t np_offset = 0;
    for (int64_t i = 0; i < node_count; ++i) {
#ifdef VERBOSE_DEBUG
        if (i % 1000 == 0) cerr << i << " of " << node_count << " ~ " << (float)i/(float)node_count * 100 << "%" << "\r";
#endif
        np_bv[np_offset] = 1; // mark node start
        //uint64_t idx = number_bool_packing::pack(i, false)+1;
        uint64_t has_steps = false;
        node_path_mm->for_values_of(i+1, [&](const std::tuple<uint64_t, uint64_t, uint64_t>& v) {
            np_iv[np_offset] = std::get<0>(v); // packed path and rev
            nr_iv[np_offset] = std::get<1>(v); // rank of step in path
            nx_iv[np_offset] = std::get<2>(v); // offset from path start on forward handle
            ++np_offset;
            has_steps = true;
        });
        if (!has_steps) ++np_offset; // skip over null entry
    }
    
    /*
     std::cerr << std::endl;
     std::cerr << "np_bv ";
     for (uint64_t i = 0; i < np_bv.size(); ++i) std:cerr << np_bv[i] << " ";
     std::cerr << std::endl;
     std::cerr << "np_iv " << np_iv << std::endl;
     std::cerr << "nr_iv " << nr_iv << std::endl;
     std::cerr << "nx_iv " << nx_iv << std::endl;
     */
    
#ifdef VERBOSE_DEBUG
    std::cerr << node_count << " of " << node_count << " ~ 100.0000%" << std::endl;
#endif
    node_path_mm.reset();
    std::remove(node_path_idx.c_str());
    // TODO, evaluate if these should be compressed more intensely
    sdsl::util::bit_compress(np_iv);
    sdsl::util::bit_compress(nr_iv);
    sdsl::util::bit_compress(nx_iv);
    sdsl::util::assign(np_bv_select, sdsl::bit_vector::select_1_type(&np_bv));
}

void XG::to_gfa(std::ostream& out) const {
    out << "H\tVN:Z:1.0" << std::endl;
    // for each node
    for_each_handle([&out,this](const handle_t& h) {
            nid_t node_id = get_id(h);
            out << "S\t" << node_id << "\t" << get_sequence(h) << std::endl;
            follow_edges(h, false, [&](const handle_t& o) {
                    out << "L\t" << node_id << "\t"
                        << (get_is_reverse(h)?"-":"+")
                        << "\t" << get_id(o) << "\t"
                        << (get_is_reverse(o)?"-":"+")
                        << "\t0M" << std::endl;
                });
            follow_edges(h, true, [&](const handle_t& o) {
                    if (get_is_reverse(o) || h == o) {
                        out << "L\t" << node_id << "\t"
                            << (get_is_reverse(h)?"-":"+")
                            << "\t" << get_id(o) << "\t"
                            << (get_is_reverse(o)?"-":"+")
                            << "\t0M" << std::endl;
                    }
                });
        });
    for_each_path_handle([&out,this](const path_handle_t& p) {
            //step_handle_t step = path_begin(p);
            out << "P\t" << get_path_name(p) << "\t";
            for_each_step_in_path(p, [this,&out](const step_handle_t& step) {
                    handle_t h = get_handle_of_step(step);
                    out << get_id(h) << (get_is_reverse(h)?"-":"+");
                    if (has_next_step(step)) out << ",";
                });
            out << "\t";
            for_each_step_in_path(p, [this,&out](const step_handle_t& step) {
                    out << get_length(get_handle_of_step(step)) << "M";
                    if (has_next_step(step)) out << ",";
                });
            out << std::endl;
        });
}

char XG::pos_char(nid_t id, bool is_rev, size_t off) const {
    assert(off < get_length(get_handle(id)));
    if (!is_rev) {
        size_t rank = id_to_rank(id);
        size_t pos = s_bv_select(rank) + off;
        assert(pos < s_iv.size());
        char c = revdna3bit(s_iv[pos]);
        return c;
    } else {
        size_t rank = id_to_rank(id);
        size_t pos = s_bv_select(rank+1) - (off+1);
        assert(pos < s_iv.size());
        char c = revdna3bit(s_iv[pos]);
        return reverse_complement(c);
    }
}

std::string XG::pos_substr(nid_t id, bool is_rev, size_t off, size_t len) const {
    if (!is_rev) {
        size_t rank = id_to_rank(id);
        size_t start = s_bv_select(rank) + off;
        assert(start < s_iv.size());
        // get until the end position, or the end of the node, which ever is first
        size_t end;
        if (!len) {
            end = s_bv_select(rank+1);
        } else {
            end = min(start + len, (size_t)s_bv_select(rank+1));
        }
        assert(end < s_iv.size());
        string s; s.resize(end-start);
        for (size_t i = start; i < s_bv.size() && i < end; ++i) {
            s[i-start] = revdna3bit(s_iv[i]);
        }
        return s;
    } else {
        size_t rank = id_to_rank(id);
        size_t end = s_bv_select(rank+1) - off;
        assert(end < s_iv.size());
        // get until the end position, or the end of the node, which ever is first
        size_t start;
        if (len > end || !len) {
            start = s_bv_select(rank);
        } else {
            start = max(end - len, (size_t)s_bv_select(rank));
        }
        assert(end < s_iv.size());
        string s; s.resize(end-start);
        for (size_t i = start; i < s_bv.size() && i < end; ++i) {
            s[i-start] = revdna3bit(s_iv[i]);
        }
        return reverse_complement(s);
    }
}

edge_t XG::edge_from_encoding(const nid_t& from, const nid_t& to, int type) const {
    bool from_rev = false;
    bool to_rev = false;
    switch (type) {
    case EDGE_TYPE_END_START:
        break;
    case EDGE_TYPE_END_END:
        to_rev = true;
        break;
    case EDGE_TYPE_START_START:
        from_rev = true;
        break;
    case EDGE_TYPE_START_END:
        from_rev = true;
        to_rev = true;
        break;
    default:
        throw std::runtime_error("Invalid edge type " + std::to_string(type) + " encountered");
        break;
    }
    return make_pair(get_handle(from, from_rev), get_handle(to, to_rev));
}

size_t XG::edge_index(const edge_t& edge) const {
    handle_t from_handle = edge.first;
    handle_t to_handle = edge.second;
    if (get_is_reverse(from_handle) && get_is_reverse(to_handle)) {
        from_handle = flip(from_handle);
        to_handle = flip(to_handle);
    }
    nid_t id = get_id(from_handle);
    size_t g = g_bv_select(id_to_rank(id));
    int edges_to_count = g_iv[g+G_NODE_TO_COUNT_OFFSET];
    int edges_from_count = g_iv[g+G_NODE_FROM_COUNT_OFFSET];
    int64_t f = g + G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * edges_to_count;
    int64_t e = f + G_EDGE_LENGTH * edges_from_count;
    int i = 1;
    for (int64_t j = f; j < e; ++i) {
        int64_t to = g+g_iv[j++];
        int type = g_iv[j++];
        edge_t curr = edge_from_encoding(id,
                                         (nid_t)g_iv[to+G_NODE_ID_OFFSET],
                                         type);
        if (curr.second == to_handle
            && get_is_reverse(curr.first) == get_is_reverse(from_handle)
            && get_is_reverse(curr.second) == get_is_reverse(to_handle)) {
            return g + i;
        }
    }
    throw std::runtime_error("Cound not find index of edge connecting " +
        std::to_string(get_id(from_handle)) + " and " + std::to_string(get_id(to_handle)));
    return 0;
}

size_t XG::get_g_iv_size(void) const {
    return g_iv.size();
}

size_t XG::id_to_rank(const nid_t& id) const {
    size_t x = id-min_id;
    if (x < 0 || x >= r_iv.size()) return 0;
    return r_iv[x];
}

size_t XG::handle_rank(const handle_t& handle) const {
    return id_to_rank(get_id(handle));
}

nid_t XG::rank_to_id(const size_t& rank) const {
    if(rank == 0) {
        cerr << "[xg] error: Request for id of rank 0" << endl;
        exit(1);
    }
    if(rank > node_count) {
        cerr << "[xg] error: Request for id of rank " << rank << "/" << node_count << endl;
        exit(1);
    }
    return g_iv[g_bv_select(rank)];
}

handle_t XG::get_handle(const nid_t& node_id, bool is_reverse) const {
    // Handles will be g vector index with is_reverse in the low bit
    
    // What rank are we looking for?
    size_t node_rank = id_to_rank(node_id);
    
    if (node_rank == 0) {
        throw runtime_error("Attempted to get handle for node " + std::to_string(node_id) + " not present in graph");
    }
    
    // Where in the g vector do we need to be
    uint64_t g = g_bv_select(node_rank);
    
    if (g + G_NODE_HEADER_LENGTH > g_iv.size()) {
        throw runtime_error("Handle for node " + std::to_string(node_id) + " with g vector offset " +
            std::to_string(g) + " is too close to end of g vector at " + std::to_string(g_iv.size()));
    }
    
    // And set the high bit if it's reverse
    return handlegraph::number_bool_packing::pack(g, is_reverse);
}

nid_t XG::get_id(const handle_t& handle) const {
    // Go get the g offset and then look up the noder ID
    return g_iv[handlegraph::number_bool_packing::unpack_number(handle) + G_NODE_ID_OFFSET];
}

bool XG::has_node(nid_t id) const {
    return id_to_rank(id) != 0;
}

bool XG::get_is_reverse(const handle_t& handle) const {
    return handlegraph::number_bool_packing::unpack_bit(handle);
}

handle_t XG::flip(const handle_t& handle) const {
    return handlegraph::number_bool_packing::toggle_bit(handle);
}

size_t XG::get_length(const handle_t& handle) const {
    return g_iv[handlegraph::number_bool_packing::unpack_number(handle) + G_NODE_LENGTH_OFFSET];
}

string XG::get_sequence(const handle_t& handle) const {
    
    // Figure out how big it should be
    size_t sequence_size = get_length(handle);
    // Allocate the sequence string
    string sequence(sequence_size, '\0');
    // Extract the node record start
    size_t g = handlegraph::number_bool_packing::unpack_number(handle);
    // Figure out where the sequence starts
    size_t sequence_start = g_iv[g + G_NODE_SEQ_START_OFFSET];
    for (int64_t i = 0; i < sequence_size; i++) {
        // Blit the sequence out
        sequence[i] = revdna3bit(s_iv[sequence_start + i]);
    }
    
    if (handlegraph::number_bool_packing::unpack_bit(handle)) {
        reverse_complement_in_place(sequence);
    }
    
    return sequence;
}

char XG::get_base(const handle_t& handle, size_t index) const {
    
    // Figure out where the sequence starts
    size_t sequence_start = g_iv[handlegraph::number_bool_packing::unpack_number(handle) + G_NODE_SEQ_START_OFFSET];
    
    // get the character
    if (get_is_reverse(handle)) {
        return reverse_complement(revdna3bit(s_iv[sequence_start + get_length(handle) - index - 1]));
    }
    else {
        return revdna3bit(s_iv[sequence_start + index]);
    }
}

string XG::get_subsequence(const handle_t& handle, size_t index, size_t size) const {
    
    // Figure out how big it should be
    size_t sequence_size = get_length(handle);
    // don't go past the end of the sequence
    size = min(size, sequence_size - index);
    // Figure out where the sequence starts
    size_t sequence_start = g_iv[handlegraph::number_bool_packing::unpack_number(handle) + G_NODE_SEQ_START_OFFSET];

    // Allocate the sequence string
    string subsequence(size, '\0');
    // unpack the sequence and handle orientation
    if (get_is_reverse(handle)) {
        for (size_t i = 0, subseq_start = sequence_start + get_length(handle) - index - size; i < size; ++i) {
            subsequence[subsequence.size() - i - 1] = reverse_complement(revdna3bit(s_iv[subseq_start + i]));
        }
    }
    else {
        for (size_t i = 0, subseq_start = sequence_start + index; i < size; ++i) {
            subsequence[i] = revdna3bit(s_iv[subseq_start + i]);
        }
    }
    return subsequence;
}

int XG::edge_type(const handle_t& from, const handle_t& to) const {
    if (get_is_reverse(from) && get_is_reverse(to)) {
        return EDGE_TYPE_START_END;
    } else if (get_is_reverse(from)) {
        return EDGE_TYPE_START_START;
    } else if (get_is_reverse(to)) {
        return EDGE_TYPE_END_END;
    } else {
        return EDGE_TYPE_END_START;
    }
}

bool XG::edge_filter(int type, bool is_to, bool want_left, bool is_reverse) const {
    // Return true if we want an edge of the given type, where we are the from
    // or to node (according to is_to), when we are looking off the right or
    // left side of the node (according to want_left), and when the node is
    // forward or reverse (accoridng to is_reverse).
    
    // First compute what we want looking off the right of a node in the forward direction.
    bool wanted = (!is_to && (type == EDGE_TYPE_END_START || type == EDGE_TYPE_END_END)) || 
        (is_to && (type == EDGE_TYPE_END_END || type == EDGE_TYPE_START_END));
    
    // We computed whether we wanted it assuming we were looking off the right. The complement is what we want looking off the left.
    wanted = wanted != want_left;
    
    // We computed whether we wanted ot assuming we were in the forward orientation. The complement is what we want in the reverse orientation.
    wanted = wanted != is_reverse;
    
    return wanted;
}

bool XG::do_edges(const size_t& g, const size_t& start, const size_t& count, bool is_to,
bool want_left, bool is_reverse, const function<bool(const handle_t&)>& iteratee) const {
    
    // OK go over all those edges
    for (size_t i = 0; i < count; i++) {
        // What edge type is the edge?
        int type = g_iv[start + i * G_EDGE_LENGTH + G_EDGE_TYPE_OFFSET];
        
        // Make sure we got a valid edge type and we haven't wandered off into non-edge data.
        assert(type >= EDGE_TYPE_MIN);
        assert(type <= EDGE_TYPE_MAX);
        
        if (edge_filter(type, is_to, want_left, is_reverse)) {
            
            // What's the offset to the other node?
            int64_t offset = g_iv[start + i * G_EDGE_LENGTH + G_EDGE_OFFSET_OFFSET];
            
            // Make sure we haven't gone off the rails into non-edge data.
            assert((int64_t) g + offset >= 0);
            assert(g + offset < g_iv.size());
            
            // Should we invert?
            // We only invert if we cross an end to end edge. Or a start to start edge
            bool new_reverse = is_reverse != (type == EDGE_TYPE_END_END || type == EDGE_TYPE_START_START);
            
            // Compose the handle for where we are going
            handle_t next_handle = handlegraph::number_bool_packing::pack(g + offset, new_reverse);
            
            // We want this edge
            
            if (!iteratee(next_handle)) {
                // Stop iterating
                return false;
            }
        }
        else {
            // TODO: delete this after using it to debug
            int64_t offset = g_iv[start + i * G_EDGE_LENGTH + G_EDGE_OFFSET_OFFSET];
            bool new_reverse = is_reverse != (type == EDGE_TYPE_END_END || type == EDGE_TYPE_START_START);
            handle_t next_handle = handlegraph::number_bool_packing::pack(g + offset, new_reverse);
        }
    }
    // Iteratee didn't stop us.
    return true;
}

bool XG::follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const {

    // Unpack the handle
    size_t g = handlegraph::number_bool_packing::unpack_number(handle);
    bool is_reverse = handlegraph::number_bool_packing::unpack_bit(handle);

    // How many edges are there of each type?
    size_t edges_to_count = g_iv[g + G_NODE_TO_COUNT_OFFSET];
    size_t edges_from_count = g_iv[g + G_NODE_FROM_COUNT_OFFSET];
    
    // Where does each edge run start?
    size_t to_start = g + G_NODE_HEADER_LENGTH;
    size_t from_start = g + G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * edges_to_count;
    
    // We will look for all the edges on the appropriate side, which means we have to check the from and to edges
    if (do_edges(g, to_start, edges_to_count, true, go_left, is_reverse, iteratee)) {
        // All the edges where we're to were accepted, so do the edges where we're from
        return do_edges(g, from_start, edges_from_count, false, go_left, is_reverse, iteratee);
    } else {
        return false;
    }
}

bool XG::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel) const {
    // This shared flag lets us bail early even when running in parallel
    bool stop_early = false;
    if (parallel) {
#pragma omp parallel
        {
#pragma omp single 
            {
                // We need to do a serial scan of the g vector because each entry is variable size.
                for (size_t g = 0; g < g_iv.size() && !stop_early;) {
                    // Make it into a handle, packing it as the node ID and using 0 for orientation
                    handle_t handle = handlegraph::number_bool_packing::pack(g, false);
                
#pragma omp task firstprivate(handle)
                    {
                        // Run the iteratee, skipping empty handles that result from discontiguous ids
                        if (get_id(handle) && !iteratee(handle)) {
                            // The iteratee is bored and wants to stop.
#pragma omp atomic write
                            stop_early = true;
                        }
                    }
                    
                    // How many edges are there of each type on this record?
                    size_t edges_to_count = g_iv[g + G_NODE_TO_COUNT_OFFSET];
                    size_t edges_from_count = g_iv[g + G_NODE_FROM_COUNT_OFFSET];
                    
                    // This record is the header plus all the edge records it contains.
                    // Decode the entry size in the same thread doing the iteration.
                    g += G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * (edges_to_count + edges_from_count);
                }
            }
            
            // The end of the single block waits for all the tasks 
        }
    } else {
        for (size_t g = 0; g < g_iv.size() && !stop_early;) {
            // Make it into a handle, packing it as the node ID and using 0 for orientation
            handle_t handle = handlegraph::number_bool_packing::pack(g, false);
            
            // Run the iteratee in-line, skipping empty handles that result from discontiguous ids
            if (get_id(handle) && !iteratee(handle)) {
                // The iteratee is bored and wants to stop.
                stop_early = true;
            }
            
            // How many edges are there of each type on this record?
            size_t edges_to_count = g_iv[g + G_NODE_TO_COUNT_OFFSET];
            size_t edges_from_count = g_iv[g + G_NODE_FROM_COUNT_OFFSET];
            
            // This record is the header plus all the edge records it contains.
            // Decode the entry size in the same thread doing the iteration.
            g += G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * (edges_to_count + edges_from_count);
        }
    }
    
    return !stop_early;
}

size_t XG::get_path_count() const {
    return paths.size();
}

bool XG::has_path(const std::string& path_name) const {
    return as_integer(get_path_handle(path_name)) != 0;
}

path_handle_t XG::get_path_handle(const std::string& path_name) const {
    // find the name in the csa
    std::string query = start_marker + path_name + end_marker;
    auto occs = locate(pn_csa, query);
    if (occs.size() > 1) {
        cerr << "error [xg]: multiple hits for " << query << endl;
        exit(1);
    }
    if(occs.size() == 0) {
        // This path does not exist. Give back 0, which can never be a real path
        // rank.
        return as_path_handle(0);
    }
    return as_path_handle(pn_bv_rank(occs[0])+1); // step past '#'
}
    
std::string XG::get_path_name(const path_handle_t& path_handle) const {
    uint64_t rank = as_integer(path_handle);
    size_t start = pn_bv_select(rank)+1; // step past '#'
    size_t end = rank == path_count ? pn_iv.size() : pn_bv_select(rank+1);
    end -= 1;  // step before '$'
    string name; name.resize(end-start);
    for (size_t i = start; i < end; ++i) {
        name[i-start] = pn_iv[i];
    }
    return name;
}

bool XG::get_is_circular(const path_handle_t& path_handle) const {
    return paths[as_integer(path_handle) - 1]->is_circular;
}

size_t XG::get_step_count(const path_handle_t& path_handle) const {
    return paths[as_integer(path_handle) - 1]->handles.size();
}
    
size_t XG::get_path_length(const path_handle_t& path_handle) const {
    return paths[as_integer(path_handle) - 1]->offsets.size();
}

handle_t XG::get_handle_of_step(const step_handle_t& step_handle) const {
    const auto& xgpath = *paths[as_integer(get_path_handle_of_step(step_handle)) - 1];
    return xgpath.handle(as_integers(step_handle)[1]);
}

path_handle_t XG::get_path_handle_of_step(const step_handle_t& step_handle) const {
    return as_path_handle(as_integers(step_handle)[0]);
}

step_handle_t XG::path_begin(const path_handle_t& path_handle) const {
    step_handle_t step;
    as_integers(step)[0] = as_integer(path_handle);
    as_integers(step)[1] = 0;
    return step;
}

step_handle_t XG::path_end(const path_handle_t& path_handle) const {
    step_handle_t step;
    as_integers(step)[0] = as_integer(path_handle);
    as_integers(step)[1] = get_step_count(path_handle);
    return step;
}

step_handle_t XG::path_back(const path_handle_t& path_handle) const {
    step_handle_t step;
    as_integers(step)[0] = as_integer(path_handle);
    as_integers(step)[1] = get_step_count(path_handle) - 1;
    return step;
    
}

step_handle_t XG::path_front_end(const path_handle_t& path_handle) const {
    step_handle_t step;
    as_integers(step)[0] = as_integer(path_handle);
    as_integers(step)[1] = -1;
    return step;
}

bool XG::has_next_step(const step_handle_t& step_handle) const {
    return (as_integers(step_handle)[1] + 1 < get_step_count(get_path_handle_of_step(step_handle))
            || (get_is_circular(get_path_handle_of_step(step_handle))
                && get_step_count(get_path_handle_of_step(step_handle)) > 0));
}

bool XG::has_previous_step(const step_handle_t& step_handle) const {
    return (as_integers(step_handle)[1] > 0
            || (get_is_circular(get_path_handle_of_step(step_handle))
                && get_step_count(get_path_handle_of_step(step_handle)) > 0));
}

step_handle_t XG::get_next_step(const step_handle_t& step_handle) const {
    step_handle_t next_step;
    as_integers(next_step)[0] = as_integers(step_handle)[0];
    as_integers(next_step)[1] = as_integers(step_handle)[1] + 1;
    if (get_is_circular(get_path_handle_of_step(step_handle))) {
        if (as_integers(next_step)[1] == get_step_count(get_path_handle_of_step(step_handle))) {
            as_integers(next_step)[1] = 0;
        }
    }
    return next_step;
}

step_handle_t XG::get_previous_step(const step_handle_t& step_handle) const {
    step_handle_t prev_step;
    as_integers(prev_step)[0] = as_integers(step_handle)[0];
    if (get_is_circular(get_path_handle_of_step(step_handle)) && as_integers(step_handle)[1] == 0) {
        as_integers(prev_step)[1] = get_step_count(get_path_handle_of_step(step_handle)) - 1;
    }
    else {
        as_integers(prev_step)[1] = as_integers(step_handle)[1] - 1;
    }
    return prev_step;

}

bool XG::for_each_path_handle_impl(const function<bool(const path_handle_t&)>& iteratee) const {
    for (size_t i = 0; i < paths.size(); i++) {
        // convert to 1-based rank
        path_handle_t path_handle = as_path_handle(i + 1);
        // execute function
        if (!iteratee(path_handle)) {
            return false;
        }
    }
    return true;
}

bool XG::for_each_step_on_handle_impl(const handle_t& handle, const function<bool(const step_handle_t&)>& iteratee) const {
    size_t off = np_bv_select(id_to_rank(get_id(handle)));
    size_t i = off;
    while (i < np_bv.size() && (off == i && np_iv[i] != 0 || np_bv[i] == 0)) {
        step_handle_t step_handle;
        as_integers(step_handle)[0] = number_bool_packing::unpack_number(as_handle(np_iv[i]));
        as_integers(step_handle)[1] = nr_iv[i];
        if (!iteratee(step_handle)) {
            return false;
        }
        ++i;
    }
    return true;
}

bool XG::for_each_step_position_on_handle(const handle_t& handle, const std::function<bool(const step_handle_t&, const bool&, const uint64_t&)>& iteratee) const {
    size_t off = np_bv_select(id_to_rank(get_id(handle)));
    size_t i = off;
    while (i < np_bv.size() && (off == i && np_iv[i] != 0 || np_bv[i] == 0)) {
        handle_t path_and_rev = as_handle(np_iv[i]);
        step_handle_t step_handle;
        as_integers(step_handle)[0] = number_bool_packing::unpack_number(path_and_rev);
        as_integers(step_handle)[1] = nr_iv[i];
        if (!iteratee(step_handle, number_bool_packing::unpack_bit(path_and_rev), nx_iv[i])) {
            return false;
        }
        ++i;
    }
    return true;
}

/// Gets the position of a given step in the path it's from
size_t XG::get_position_of_step(const step_handle_t& step) const {
    const auto& xgpath = *paths[as_integer(get_path_handle_of_step(step)) - 1];
    return xgpath.handle_start(as_integers(step)[1]);
}

/// Get the step at a given position
step_handle_t XG::get_step_at_position(const path_handle_t& path, const size_t& position) const {
    
    if (position >= get_path_length(path)) {
        throw runtime_error("Cannot get position " + std::to_string(position) + " along path " +
            get_path_name(path) + " of length " + std::to_string(get_path_length(path)));
    }
    
    const auto& xgpath = *paths[as_integer(path) - 1];
    step_handle_t step;
    as_integers(step)[0] = as_integer(path);
    as_integers(step)[1] = xgpath.step_rank_at_position(position);
    return step;
}

size_t XG::get_node_count() const {
    return this->node_count;
}

nid_t XG::min_node_id() const {
    return min_id;
}
    
nid_t XG::max_node_id() const {
    return max_id;
}

size_t XG::max_node_rank(void) const {
    return s_bv_rank(s_bv.size());
}

nid_t XG::node_at_vector_offset(const size_t& offset) const {
    return rank_to_id(s_bv_rank(offset));
}

size_t XG::node_vector_offset(const nid_t& id) const {
    return s_bv_select(id_to_rank(id));
}

size_t XG::node_graph_idx(const nid_t& id) const {
    return g_bv_select(id_to_rank(id));
}

bool XG::path_contains_handle(const std::string& name, const handle_t& handle) const {
    bool contains_handle = false;
    for_each_step_position_on_handle(handle, [&](const step_handle_t& step, const bool& is_rev, const uint64_t& pos) {
            path_handle_t path = as_path_handle(as_integers(step)[0]);
            if (get_path_name(path) == name) {
                contains_handle = true;
                return false; // break
            } else {
                return true; // keep going
            }
        });
    return false;
}

std::vector<path_handle_t> XG::paths_of_handle(const handle_t& handle) const {
    std::vector<path_handle_t> path_handles;
    for_each_step_position_on_handle(handle, [&](const step_handle_t& step, const bool& is_rev, const uint64_t& pos) {
            path_handles.push_back(as_path_handle(as_integers(step)[0]));
            return true;
        });
    return path_handles;
}

std::pair<pos_t, int64_t> XG::next_path_position(const pos_t& pos, const int64_t& max_search) const {
    handle_t h_fwd = get_handle(id(pos), is_rev(pos));
    handle_t h_rev = get_handle(id(pos), !is_rev(pos));
    int64_t fwd_seen = offset(pos);
    int64_t rev_seen = get_length(h_fwd) - offset(pos);
    std::pair<pos_t, int64_t> fwd_next = make_pair(make_pos_t(0,false,0), std::numeric_limits<int64_t>::max());
    std::pair<pos_t, int64_t> rev_next = make_pair(make_pos_t(0,false,0), std::numeric_limits<int64_t>::max());
    follow_edges(h_fwd, false, [&](const handle_t& n) {
            nid_t id = get_id(n);
            if (!paths_of_handle(n).empty()) {
                fwd_next = make_pair(make_pos_t(id, get_is_reverse(n), 0), fwd_seen);
                return false;
            } else {
                fwd_seen += get_length(n);
                return fwd_seen < max_search;
            }
        });
    follow_edges(h_rev, false, [&](const handle_t& n) {
            nid_t id = get_id(n);
            if (!paths_of_handle(n).empty()) {
                rev_next = make_pair(make_pos_t(id, !get_is_reverse(n), 0), rev_seen);
                return false;
            } else {
                rev_seen += get_length(n);
                return rev_seen < max_search;
            }
        });
    if (fwd_next.second <= rev_next.second) {
        return fwd_next;
    } else {
        rev_next.second = -rev_next.second;
        return rev_next;
    }
}

std::pair<nid_t, std::vector<path_handle_t> > XG::nearest_path_node(const nid_t& id, int64_t max_steps) const {
    std::set<nid_t> todo;
    std::set<nid_t> seen;
    todo.insert(id);
    int i = 0;
    while (!todo.empty() && i++ < max_steps) {
        std::set<nid_t> next;
        for (auto& id : todo) {
            if (seen.count(id)) continue;
            seen.insert(id);
            std::vector<path_handle_t> path_ids = paths_of_handle(get_handle(id, false));
            if (!path_ids.empty()) {
                // if we found a node on a path, return
                return make_pair(id, path_ids);
            } else {
                // collect neighboring nodes
                handle_t handle = get_handle(id, false);
                follow_edges(handle, false, [&](const handle_t& n) {
                        next.insert(get_id(n));
                    });
                follow_edges(handle, true, [&](const handle_t& p) {
                        next.insert(get_id(p));
                    });
            }
        }
        todo = next;
    }
    return make_pair(id, std::vector<path_handle_t>());
}

// like above, but find minumum over list of paths.  if names is empty, do all paths
// don't actually take strict minumum over all paths.  rather, prefer paths that
// contain the nodes when possible. 
int64_t XG::min_approx_path_distance(const nid_t& id1, const nid_t& id2) const {
    int64_t min_distance = numeric_limits<int64_t>::max();
    std::pair<nid_t, std::vector<path_handle_t> > near1 = nearest_path_node(id1);
    std::pair<nid_t, std::vector<path_handle_t> > near2 = nearest_path_node(id2);
    if (near1.second.size() || near2.second.size()) {
        std::unordered_map<path_handle_t, size_t> paths;
        for (auto& i : near1.second) paths[i]++;
        for (auto& i : near2.second) paths[i]++;
        for (auto& i : paths) {
            if (i.second < 2) continue;
            //auto name = get_path_name(i.first);
            auto pos1s = position_in_path(get_handle(near1.first, false), i.first);
            auto pos2s = position_in_path(get_handle(near2.first, false), i.first);
            for (auto& p1 : pos1s) {
                for (auto& p2 : pos2s) {
                    int64_t distance = abs((int64_t)p1 - (int64_t)p2);
                    min_distance = min(distance, min_distance);
                }
            }
        }
    }
    return min_distance;
}

void XG::for_path_range(const std::string& name, int64_t start, int64_t stop,
                        std::function<void(const handle_t&)> lambda, bool is_rev) const {
    // what is the node at the start, and at the end
    auto& path = *paths[as_integer(get_path_handle(name))-1];
    size_t plen = path.offsets.size();
    if (start > plen) return; // no overlap with path
    // careful not to exceed the path length
    if (stop >= plen) stop = plen-1;
    if (is_rev) {
        start = plen - start;
        stop = plen - stop;
    }
    size_t pr1 = path.offsets_rank(start+1)-1;
    size_t pr2 = path.offsets_rank(stop+1)-1;
    // Grab the IDs visited in order along the path
    for (size_t i = pr1; i <= pr2; ++i) {
        lambda(path.handle(i));
    }
}

std::vector<size_t> XG::position_in_path(const handle_t& handle, const path_handle_t& path) const {
    std::vector<size_t> pos_in_path;
    size_t path_length = get_path_length(path);
    for_each_step_position_on_handle(handle, [&](const step_handle_t& step, const bool& is_rev, const uint64_t& pos) {
            path_handle_t path_handle = as_path_handle(as_integers(step)[0]);
            if (path_handle == path) {
                pos_in_path.push_back(get_is_reverse(handle) ? path_length - pos - get_length(handle) : pos);
            }
            return true;
        });
    return pos_in_path;
}

std::unordered_map<path_handle_t, std::vector<size_t> > XG::position_in_paths(const handle_t& handle, const size_t& offset) const {
    std::unordered_map<path_handle_t, std::vector<size_t> > positions;
    for_each_step_position_on_handle(handle, [&](const step_handle_t& step, const bool& is_rev, const uint64_t& pos) {
            path_handle_t path = as_path_handle(as_integers(step)[0]);
            size_t path_length = get_path_length(path);
            positions[path].push_back(get_is_reverse(handle) ? path_length - pos - get_length(handle) : pos);
            return true;
        });
    return positions;
}

std::unordered_map<path_handle_t, std::vector<std::pair<size_t, bool> > > XG::offsets_in_paths(const pos_t& gpos) const {
    std::unordered_map<path_handle_t, std::vector<std::pair<size_t, bool> > > positions;
    handle_t handle = get_handle(id(gpos), is_rev(gpos));
    size_t handle_length = get_length(handle);
    for_each_step_position_on_handle(handle, [&](const step_handle_t& step, const bool& rev, const uint64_t& pos) {
            path_handle_t path = as_path_handle(as_integers(step)[0]);
            size_t path_length = get_path_length(path);
            bool dir = rev != is_rev(gpos);
            size_t node_forward_strand_offset = is_rev(gpos) ? (handle_length - offset(gpos) - 1) : offset(gpos);
            size_t off = pos + (rev ?
                                (handle_length - node_forward_strand_offset - 1) :
                                node_forward_strand_offset);
            positions[path].push_back(std::make_pair(off, dir));
            return true;
        });
    return positions;
}

std::unordered_map<path_handle_t, std::vector<std::pair<size_t, bool> > > XG::nearest_offsets_in_paths(const pos_t& pos, int64_t max_search) const {
    std::pair<pos_t, int64_t> pz = next_path_position(pos, max_search);
    auto& path_pos = pz.first;
    auto& diff = pz.second;
    if (id(path_pos)) {
        // TODO apply approximate offset, second in pair returned by next_path_position
        auto offsets = offsets_in_paths(path_pos);
        for (auto& o : offsets) {
            for (auto& p : o.second) {
                p.first += diff;
            }
        }
        return offsets;
    } else {
        std::unordered_map<path_handle_t, std::vector<std::pair<size_t, bool> > > empty;
        return empty;
    }
}

handle_t XG::handle_at_path_position(const path_handle_t& path, size_t pos) const {
    return paths[as_integer(path)]->handle_at_position(pos);
}

size_t XG::node_start_at_path_position(const path_handle_t& path, size_t pos) const {
    size_t p = as_integer(path);
    size_t position_rank = paths[p]->offsets_rank(pos+1);
    return paths[p]->offsets_select(position_rank);
}

pos_t XG::graph_pos_at_path_position(const path_handle_t& path_handle, size_t path_pos) const {
    auto& path = *paths[as_integer(path_handle)];
    //handle_t path.handle_at_position(path_pos);
    // what's the path offset in the path?
    size_t x = path.step_rank_at_position(path_pos);
    handle_t handle = path.handle(x);
    return make_pos_t(get_id(handle), get_is_reverse(handle), path_pos - path.handle_start(x));
}

namespace temp_file {

// We use this to make the API thread-safe
std::recursive_mutex monitor;

std::string temp_dir;

/// Because the names are in a static object, we can delete them when
/// std::exit() is called.
struct Handler {
    std::set<std::string> filenames;
    std::string parent_directory;
    ~Handler() {
        // No need to lock in static destructor
        for (auto& filename : filenames) {
            std::remove(filename.c_str());
        }
        if (!parent_directory.empty()) {
            // There may be extraneous files in the directory still (like .fai files)
            auto directory = opendir(parent_directory.c_str());
            
            dirent* dp;
            while ((dp = readdir(directory)) != nullptr) {
                // For every item still in it, delete it.
                // TODO: Maybe eventually recursively delete?
                std::remove((parent_directory + "/" + dp->d_name).c_str());
            }
            closedir(directory);
            
            // Delete the directory itself
            std::remove(parent_directory.c_str());
        }
    }
} handler;

std::string create(const std::string& base) {
    std::lock_guard<recursive_mutex> lock(monitor);

    if (handler.parent_directory.empty()) {
        // Make a parent directory for our temp files
        string tmpdirname_cpp = get_dir() + "/xg-XXXXXX";
        char* tmpdirname = new char[tmpdirname_cpp.length() + 1];
        strcpy(tmpdirname, tmpdirname_cpp.c_str());
        auto got = mkdtemp(tmpdirname);
        if (got != nullptr) {
            // Save the directory we got
            handler.parent_directory = got;
        } else {
            cerr << "[xg]: couldn't create temp directory: " << tmpdirname << endl;
            exit(1);
        }
        delete [] tmpdirname;
    }

    std::string tmpname = handler.parent_directory + "/" + base + "XXXXXX";
    // hack to use mkstemp to get us a safe temporary file name
    int fd = mkstemp(&tmpname[0]);
    if(fd != -1) {
        // we don't leave it open; we are assumed to open it again externally
        close(fd);
    } else {
        cerr << "[xg]: couldn't create temp file on base "
             << base << " : " << tmpname << endl;
        exit(1);
    }
    handler.filenames.insert(tmpname);
    return tmpname;
}

std::string create() {
    // No need to lock as we call this thing that locks
    return create("xg-");
}

void remove(const std::string& filename) {
    std::lock_guard<recursive_mutex> lock(monitor);
    
    std::remove(filename.c_str());
    handler.filenames.erase(filename);
}

void set_dir(const std::string& new_temp_dir) {
    std::lock_guard<recursive_mutex> lock(monitor);
    
    temp_dir = new_temp_dir;
}

std::string get_dir() {
    std::lock_guard<recursive_mutex> lock(monitor);

    // Get the default temp dir from environment variables.
    if (temp_dir.empty()) {
        const char* system_temp_dir = nullptr;
        for(const char* var_name : {"TMPDIR", "TMP", "TEMP", "TEMPDIR", "USERPROFILE"}) {
            if (system_temp_dir == nullptr) {
                system_temp_dir = getenv(var_name);
            }
        }
        temp_dir = (system_temp_dir == nullptr ? "/tmp" : system_temp_dir);
    }

    return temp_dir;
}

}

}
