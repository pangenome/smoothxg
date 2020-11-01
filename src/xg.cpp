#include "xg.hpp"

//#include "ips4o.hpp"
#include "mmmultimap.hpp"

#include <bitset>
#include <arpa/inet.h>
#include <mutex>

#include <handlegraph/util.hpp>

#include "gfakluge.hpp"

//#define VERBOSE_DEBUG
//#define debug_path_index
//#define DEBUG_CONSTRUCTION
//#define debug_print_graph

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

uint32_t XG::get_magic_number(void) const {
    // This is not a nice string
    return 4143290017ul;
}
    
void XG::deserialize_members(std::istream& in) {

    if (!in.good()) {
        throw XGFormatError("Index file does not exist or index stream cannot be read");
    }

    // Version 0 is the last XG format without an explicit version specifier.
    // If we find a version specifier we will up this.
    uint32_t file_version = 0;
    
    // Magic value handling is now external.

    // We need to look for the old magic value (v13 or lower)
    bool have_old_magic = false;
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
            
            // Remember we had the old magic value
            have_old_magic = true;
            
        } else {
            // Put back both characters
            in.unget();
            in.unget();
        }        
    } else {
        // Put back the one character
        in.unget();
    }
    
    if (!have_old_magic) {
        // New SerializableHandleGraph magic is done for us
        
        // First 4 bytes we see are version number
        in.read((char*) &file_version, sizeof(file_version));
        // Make sure to convert from network to host byte order
        file_version = ntohl(file_version);
    }
    
    if (have_old_magic && file_version > 13) {
        // This shouldn't happen. If it does happen, we must be reading
        // something wrong (new style version number started with X and G
        // bytes)?
        throw XGFormatError("XG index file version " + std::to_string(file_version) +
                            " has old-style magic number.");
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
        case 13:
        case 14:
        case 15:
            std::cerr << "warning:[XG] Loading an out-of-date XG format. "
                      << "For better performance over repeated loads, consider recreating this XG index." << std::endl;
            // Fall through
        case 16:
            {
                sdsl::read_member(seq_length, in);
                sdsl::read_member(node_count, in);
                sdsl::read_member(edge_count, in);
                sdsl::read_member(path_count, in);
                size_t entity_count = node_count + edge_count;
                sdsl::read_member(min_id, in);
                sdsl::read_member(max_id, in);
                
                if (file_version <= 8) {
                    // Load the old id int vector to skip
                    sdsl::int_vector<> i_iv;
                    i_iv.load(in);
                }
                r_iv.load(in);
                
                // if we've rejiggered the offsets in the g vector we need to
                // hold onto the rank vector to translate offsets throughout
                // the deserialization
                sdsl::rank_support_v<1> old_g_bv_rank;
                sdsl::bit_vector old_g_bv;
                
                if (file_version <= 14) {
                    // we need to reencode the old g vector with the new edges
                    sdsl::int_vector<> old_giv;
                    old_giv.load(in);
                    old_g_bv.load(in);
                    old_g_bv_rank.load(in, &old_g_bv);
                    {
                        // we don't actually need the select vector, load and discard it
                        sdsl::bit_vector::select_1_type old_g_bv_select;
                        old_g_bv_select.load(in, &old_g_bv);
                    }
                    reencode_old_g_vector(old_giv, old_g_bv_rank);
                }
                else {
                    // we can load the up-to-date g vector encoding
                    g_iv.load(in);
                    g_bv.load(in);
                    g_bv_rank.load(in, &g_bv);
                    g_bv_select.load(in, &g_bv);
                }

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
                    if (file_version > 12 && file_version <= 14) {
                        // the paths consist of handles, but we've changed
                        // the offsets of items in the g vector, so we need
                        // to resync
                        path->sync_offsets(old_g_bv_rank, g_bv_select);
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
#ifdef debug_print_graph
    cerr << "printing deserialized graph" << endl;
    print_graph();
#endif
}

void XGPath::load(std::istream& in) {
    sdsl::read_member(min_handle, in);
    handles.load(in);
    offsets.load(in);
    offsets_rank.load(in, &offsets);
    offsets_select.load(in, &offsets);
    sdsl::read_member(is_circular, in);    
}

void XGPath::sync_offsets(const sdsl::rank_support_v<1>& old_g_bv_rank,
                          const sdsl::bit_vector::select_1_type& g_bv_select) {
    
    // make a temporary vector to hold uncompressed handles
    sdsl::int_vector<> handles_iv;
    sdsl::util::assign(handles_iv, sdsl::int_vector<>(handles.size()));
    
    // sync the offsets and compute the minimum handle
    uint64_t new_min_handle = numeric_limits<uint64_t>::max();
    for (size_t i = 0; i < handles.size(); ++i) {
        handle_t old_handle = handle(i);
        size_t old_offset = number_bool_packing::unpack_number(old_handle);
        size_t new_offset = g_bv_select(old_g_bv_rank(old_offset) + 1);
        handle_t new_handle = number_bool_packing::pack(new_offset,
                                                        number_bool_packing::unpack_bit(old_handle));
        handles_iv[i] = as_integer(new_handle);
        new_min_handle = min<uint64_t>(new_min_handle, as_integer(new_handle));
    }
    // apply the min handle offset
    min_handle = as_handle(new_min_handle);
    for (size_t i = 0; i < handles_iv.size(); ++i) {
        handles_iv[i] -= new_min_handle;
    }
    
    // compress the new handle vector and replace the old one
    sdsl::util::assign(handles, sdsl::enc_vector<>(handles_iv));
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
            nid_t id_offset = 0;
            if (file_version >= 8) {
                // IDs are stored relative to a minimum offset
                sdsl::read_member(id_offset, in);
            }
            
            sdsl::wt_gmr<> ids;
            ids.load(in);
            
            sdsl::sd_vector<> directions;
            directions.load(in);
            // compute the minimum handle
            min_handle = handlegraph::as_handle(numeric_limits<uint64_t>::max());
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
    
    // translate the offset from a straight bit_vector to an rrr_vector
    {
        sdsl::bit_vector offsets_bv;
        offsets_bv.load(in);
        
        // skip its rank and select support
        {
            sdsl::rank_support_v<1> offsets_bv_rank;
            offsets_bv_rank.load(in, &offsets_bv);
        }
        {
            sdsl::bit_vector::select_1_type offsets_bv_select;
            offsets_bv_select.load(in, &offsets_bv);
        }
        
        // reencode it as the rrr_vector we want
        sdsl::util::assign(offsets, sdsl::rrr_vector<>(offsets_bv));
    }
    
    // recreate the rank and select support
    sdsl::util::assign(offsets_rank, sdsl::rrr_vector<>::rank_1_type(&offsets));
    sdsl::util::assign(offsets_select, sdsl::rrr_vector<>::select_1_type(&offsets));
    
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

#ifdef debug_path_index
    std::cerr << "Constructing xgpath for path with handles:" << std::endl;
    for (handle_t visiting : path) {
        std::cerr << "\t" << as_integer(visiting) << std::endl;
    }
#endif

    // The circularity flag is just a normal bool
    this->is_circular = is_circular;

    // handle integer values, the literal path
    sdsl::int_vector<> handles_iv;
    sdsl::util::assign(handles_iv, sdsl::int_vector<>(path.size()));
    // directions of traversal (typically forward, but we allow backwards)
    sdsl::bit_vector directions_bv;
    sdsl::util::assign(directions_bv, sdsl::bit_vector(path.size()));

    size_t path_off = 0;
    size_t members_off = 0;
    size_t positions_off = 0;
    size_t path_length = 0;
    // Start out with the max integer, as a handle, as our minimum-valued handle in the path.
    uint64_t min_handle_int = (path.size() ? as_integer(path[0]) : 0);

    // determine min handle value which occurs
    for (size_t i = 1; i < path.size(); ++i) {
        if (as_integer(path[i]) < min_handle_int) {
            min_handle_int = as_integer(path[i]);
        }
    }
    min_handle = as_handle(min_handle_int);

#ifdef debug_path_index
    std::cerr << "Basing on minimum handle value " << as_integer(min_handle) << " (aka " << min_handle_int << ")" << std::endl;
#endif
    
    // determine total length and record handles
    for (size_t i = 0; i < path.size(); ++i) {
        const handle_t& handle = path[i];
        path_length += graph.get_length(handle);
        handles_iv[i] = as_integer(local_handle(handle));
        
#ifdef debug_path_index
        std::cerr << "Recorded handle as " << handles_iv[i] << std::endl;
#endif
        
        // we will explode if the node isn't in the graph
    }
    sdsl::util::assign(handles, sdsl::enc_vector<>(handles_iv));
    
#ifdef debug_path_index
    for (size_t i = 0; i < path.size(); i++) {
        std::cerr << "Encoded handle as " << handles[i] << std::endl;
    }
#endif

    // make the bitvector for path offsets
    sdsl::bit_vector offsets_bv;
    sdsl::util::assign(offsets_bv, sdsl::bit_vector(path_length));

    for (size_t i = 0; i < path.size(); ++i) {
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
        throw std::runtime_error("Handle with value " + std::to_string(as_integer(handle)) +
            " cannot be converted to local space based at min handle with value " + std::to_string(as_integer(min_handle)));
    } else {
        return as_handle(as_integer(handle)-as_integer(min_handle));
    }
}

handle_t XGPath::external_handle(const handle_t& handle) const {
    return as_handle(as_integer(handle)+as_integer(min_handle));
}

size_t XG::serialize_and_measure(ostream& out, sdsl::structure_tree_node* s, std::string name) const {
    // TODO: this depends on SerializableHandleGraph's internals.
    uint32_t magic_number = htonl(get_magic_number());
    out.write((char*) &magic_number, sizeof(magic_number) / sizeof(char));
    return serialize_members_and_measure(out, s, name);
}

void XG::serialize_members(ostream& out) const {
    serialize_members_and_measure(out);
}

size_t XG::serialize_members_and_measure(ostream& out, sdsl::structure_tree_node* s, std::string name) const {

    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
    size_t written = 0;
    
    // Start with the version number; SerializableHandleGraph handles the magic
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
    from_enumerators(for_each_sequence, for_each_edge, for_each_path_element, validate, basename);
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
    min_id = std::numeric_limits<nid_t>::max();
    max_id = 0;
    // get information about graph size and id ranges
#ifdef VERBOSE_DEBUG
    cerr << "computing graph sequence length and node count" << endl;
#endif
    for_each_sequence([&](const std::string& seq, const nid_t& id) {
        // min id starts at 0
        min_id = std::min(min_id, id);
        max_id = std::max(max_id, id);
        if (seq.empty()) {
            // XG can't store empty nodes, because of its one-1-per-node combinde sequence lookup bit vector.
            // If the 1s collide, bad things happen.
            // Print a useful report for a user (instead of a crash for an uncaught exception) and abort everything.
            std::cerr << "error[XG::from_enumerators]: Graph contains empty node " << id << " which cannot be stored in an XG" << std::endl;
            exit(1);
        }
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
              << path_count << " paths" << std::endl
              << "node ids run from " << min_id << " to " << max_id << std::endl;
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
    std::string edge_left_side_idx = basename + ".left.mm";
    std::string edge_right_side_idx = basename + ".right.mm";
    auto edge_left_side_mm = std::make_unique<mmmulti::map<uint64_t, uint64_t>>(edge_left_side_idx);
    auto edge_right_side_mm = std::make_unique<mmmulti::map<uint64_t, uint64_t>>(edge_right_side_idx);
    edge_left_side_mm->open_writer();
    edge_right_side_mm->open_writer();
    size_t num_reversing_self_edges = 0;
    for_each_edge([&](const nid_t& from_id, const bool& from_rev, const nid_t& to_id, const bool& to_rev) {
        if (from_rev) {
            edge_left_side_mm->append(as_integer(temp_get_handle(from_id, false)),
                                      as_integer(temp_get_handle(to_id, !to_rev)));
        }
        else {
            edge_right_side_mm->append(as_integer(temp_get_handle(from_id, false)),
                                       as_integer(temp_get_handle(to_id, to_rev)));
        }
        if (from_id != to_id || from_rev == to_rev) {
            // we're not double-counting a reversing self loop
            if (!to_rev) {
                edge_left_side_mm->append(as_integer(temp_get_handle(to_id, false)),
                                          as_integer(temp_get_handle(from_id, from_rev)));
            }
            else {
                edge_right_side_mm->append(as_integer(temp_get_handle(to_id, false)),
                                           as_integer(temp_get_handle(from_id, !from_rev)));
            }
        }
        else {
            ++num_reversing_self_edges;
        }
    });
    handle_t max_handle = number_bool_packing::pack(r_iv.size(), true);
    edge_left_side_mm->index(get_thread_count(), as_integer(max_handle));
    edge_right_side_mm->index(get_thread_count(), as_integer(max_handle));

    // calculate g_iv size (header + edges stored twice, except reversing self edges)
    size_t g_iv_size = node_count * G_NODE_HEADER_LENGTH + (edge_count * 2 - num_reversing_self_edges);
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
        g_bv[g] = 1; // mark record start for later query
        g_iv[g++] = id;
        g_iv[g++] = node_vector_offset(id);
        g_iv[g++] = temp_node_size(id);
        size_t left_edge_count = 0;
        size_t right_edge_count = 0;
        size_t left_edge_count_idx = g++;
        size_t right_edge_count_idx = g++;
        handle_t edge_head = temp_get_handle(id, false);
        edge_left_side_mm->for_unique_values_of(as_integer(edge_head), [&](const uint64_t& edge_tail) {
            g_iv[g++] = edge_tail;
            ++left_edge_count;
        });
        edge_right_side_mm->for_unique_values_of(as_integer(edge_head), [&](const uint64_t& edge_tail) {
            g_iv[g++] = edge_tail;
            ++right_edge_count;
        });
        g_iv[left_edge_count_idx] = left_edge_count;
        g_iv[right_edge_count_idx] = right_edge_count;
    }

#ifdef VERBOSE_DEBUG
    std::cerr << node_count << " of " << node_count << " ~ 100.0000%" << std::endl;
#endif

    // cleanup our mmmultimap
    edge_left_side_mm.reset();
    edge_right_side_mm.reset();
    std::remove(edge_left_side_idx.c_str());
    std::remove(edge_right_side_idx.c_str());

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
        size_t edge_count = g_iv[g+G_NODE_LEFT_COUNT_OFFSET] + g_iv[g+G_NODE_RIGHT_COUNT_OFFSET];
        size_t begin = g + G_NODE_HEADER_LENGTH;
        for (int64_t j = begin, end = begin + edge_count; j < end; ++j) {
            uint64_t edge_tail = g_iv[j];
            bool tail_rev = edge_tail & 1;
            size_t tail_offset = g_bv_select(edge_tail >> 1);
            g_iv[j] = encode_edge(g, tail_offset, tail_rev);
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

#ifdef debug_path_index
        std::cerr << "Creating XGPath for path " << curr_path_name << " with visits:" << std::endl;
        for (handle_t visiting : curr_path_steps) {
             std::cerr << "\tHandle for " << get_id(visiting) << (get_is_reverse(visiting) ? "-" : "+") << " with value "
                    << as_integer(visiting) << std::endl;
        }
#endif
        size_t unique_member_count = 0;
        // check that curr_path_name doesn't have our delimiter
        // if it does, make a warning and replace with _
        if (curr_path_name.find(path_name_csa_delim) != std::string::npos) {
            std::cerr << "[xg::XG] warning, path name " << curr_path_name
                      << " contains delimiter used by xg's CSA, changing to ";
            std::replace(curr_path_name.begin(),
                         curr_path_name.end(),
                         path_name_csa_delim, '_');
            std::cerr << curr_path_name << std::endl;
        }
        path_names += path_name_csa_delim + curr_path_name;
        XGPath* path = new XGPath(curr_path_name, curr_path_steps,
                                  curr_is_circular,
                                  *this);
        paths.push_back(path);
        
#ifdef debug_path_index
        std::cerr << "paths[" << paths.size() - 1 << "] = " << curr_path_name << " @ " << path << std::endl;
#endif
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
    path_names += path_name_csa_delim; // final delimiter for search symmetry
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
        if (path_names[i] == path_name_csa_delim) {
            pn_bv[i] = 1; // register name start
        }
    }
    sdsl::util::assign(pn_bv_rank, sdsl::rank_support_v<1>(&pn_bv));
    sdsl::util::assign(pn_bv_select, sdsl::bit_vector::select_1_type(&pn_bv));
    
    // By default, SDSL uses the working directory for temporary files. Getting around it is
    // somewhat complicated.
    sdsl::cache_config config;
    config.dir = temp_file::get_dir();
    {
        sdsl::int_vector_buffer<8> text(sdsl::cache_file_name(sdsl::conf::KEY_TEXT, config), std::ios::out);
        for (char c : path_names) {
            text.push_back(c);
        }
        text.push_back(0); // CSA construction needs an endmarker.
    }
    sdsl::register_cache_file(sdsl::conf::KEY_TEXT, config);
    sdsl::construct(pn_csa, sdsl::cache_file_name(sdsl::conf::KEY_TEXT, config), config, 1);

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
//        path_ids_mb_size += size_in_mega_bytes(path->ids);
//        path_dir_mb_size += size_in_mega_bytes(path->directions);
//        path_pos_mb_size += size_in_mega_bytes(path->positions);
//        path_ranks_mb_size += size_in_mega_bytes(path->ranks);
        path_offsets_mb_size += size_in_mega_bytes(path->offsets);
    }
    cerr << "path ids size " << path_ids_mb_size << endl;
    cerr << "path directions size " << path_dir_mb_size << endl;
    cerr << "path positions size " << path_pos_mb_size << endl;
    cerr << "path ranks size " << path_ranks_mb_size << endl;
    cerr << "path offsets size " << path_offsets_mb_size << endl;

#endif

#ifdef debug_print_graph
    print_graph();
#endif

}

void XG::print_graph() const {
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
        int edges_left_count = g_iv[g+G_NODE_LEFT_COUNT_OFFSET];
        int edges_right_count = g_iv[g+G_NODE_RIGHT_COUNT_OFFSET];
        int sequence_size = g_iv[g+G_NODE_LENGTH_OFFSET];
        size_t seq_start = g_iv[g+G_NODE_SEQ_START_OFFSET];
        cerr << id << " ";
        for (int64_t j = seq_start; j < seq_start+sequence_size; ++j) {
            cerr << revdna3bit(s_iv[j]);
        } cerr << " : ";
        int64_t l = g + G_NODE_HEADER_LENGTH;
        int64_t r = g + G_NODE_HEADER_LENGTH + edges_left_count;
        cerr << " left ";
        for (int64_t j = l; j < r; ++j) {
            cerr << rank_to_id(g_bv_rank(g+edge_relative_offset(g_iv[j]))+1) << " ";
        }
        cerr << " right ";
        for (int64_t j = r; j < r + edges_right_count; ++j) {
            cerr << rank_to_id(g_bv_rank(g+edge_relative_offset(g_iv[j]))+1) << " ";
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

bool XG::edge_orientation(const uint64_t& raw_edge) {
    return raw_edge & 1;
}

int64_t XG::edge_relative_offset(const uint64_t& raw_edge) {
    return raw_edge & 2 ? -int64_t(raw_edge >> 2) - 1 : raw_edge >> 2;
}

uint64_t XG::encode_edge(const size_t& offest_from, const size_t& offset_to,
                         const bool& to_rev) {
    // first bit encodes reverse, the remaining bits encode a signed difference in
    // offsets in a small unsigned integer by alternating positive and negative values
    uint64_t edge;
    if (offset_to >= offest_from) {
        edge = (offset_to - offest_from) << 2;
    }
    else {
        edge = ((offest_from - offset_to) << 2) - 2;
    }
    edge |= (uint64_t) to_rev;
    return edge;
}

void XG::index_node_to_path(const std::string& basename) {
    
    // node -> paths
    // use the mmmultimap...
    std::string node_path_idx = basename + ".node_path.mm";
    auto node_path_mm = std::make_unique<mmmulti::map<uint64_t, std::tuple<uint64_t, uint64_t, uint64_t>>>(node_path_idx);
    node_path_mm->open_writer();
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
    node_path_mm->index(get_thread_count(), node_count+1);
    
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
            out << "S\t" << node_id << "\t" << get_sequence(h);
            std::vector<std::string> path_pos;
            for_each_step_position_on_handle(h, [&](const step_handle_t& step, const bool& rev, const uint64_t& pos) {
                    stringstream ss;
                    ss << "[\"" << get_path_name(as_path_handle(as_integers(step)[0])) << "\"," << (rev ? "\"-\"" : "\"+\"") << "," << pos << "]";
                    path_pos.push_back(ss.str());
                    return true;
                });
            if (!path_pos.empty()) {
                out << "\t" << "PO:J:[";
                out << path_pos[0];
                for (uint64_t i = 1; i < path_pos.size(); ++i) {
                    out << "," << path_pos[i];
                }
                out << "]";
            }
            out << std::endl;
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

void XG::orientation_from_old_edge_type(int type, bool& from_rev, bool& to_rev) const {
    switch (type) {
        case OLD_EDGE_TYPE_END_START:
            from_rev = false;
            to_rev = false;
            break;
        case OLD_EDGE_TYPE_END_END:
            from_rev = false;
            to_rev = true;
            break;
        case OLD_EDGE_TYPE_START_START:
            from_rev = true;
            to_rev = false;
            break;
        case OLD_EDGE_TYPE_START_END:
            from_rev = true;
            to_rev = true;
            break;
        default:
            throw std::runtime_error("Invalid edge type " + std::to_string(type) + " encountered");
            break;
    }
}

void XG::reencode_old_g_vector(const sdsl::int_vector<>& old_g_iv, const sdsl::rank_support_v<1>& old_g_bv_rank) {
    
    size_t total_num_edges = (old_g_iv.size() - G_NODE_HEADER_LENGTH * node_count) / OLD_G_EDGE_LENGTH;
    
    // compute the new size and allocate the
    size_t g_iv_size = G_NODE_HEADER_LENGTH * node_count + total_num_edges;
    sdsl::util::assign(g_iv, sdsl::int_vector<>(g_iv_size));
    sdsl::util::assign(g_bv, sdsl::bit_vector(g_iv_size));
    
#ifdef VERBOSE_DEBUG
    cerr << "converting g vector for graph with " << node_count << " nodes and " << total_num_edges << " edges" << endl;
#endif
    
    for (size_t old_g_idx = 0, new_g_idx = 0; new_g_idx < g_iv.size(); ) {
#ifdef VERBOSE_DEBUG
        cerr << "transfering node info for node " << old_g_iv[old_g_idx + G_NODE_ID_OFFSET] << " from old index " << old_g_idx << " to new index " << new_g_idx << endl;
#endif
        
        // record the new start position in the g vector
        g_bv[new_g_idx] = 1;
        // copy over the fields that are unchanged between old and new
        g_iv[new_g_idx + G_NODE_ID_OFFSET] = old_g_iv[old_g_idx + G_NODE_ID_OFFSET];
        g_iv[new_g_idx + G_NODE_SEQ_START_OFFSET] = old_g_iv[old_g_idx + G_NODE_SEQ_START_OFFSET];
        g_iv[new_g_idx + G_NODE_LENGTH_OFFSET] = old_g_iv[old_g_idx + G_NODE_LENGTH_OFFSET];
        // ignore the left/right edge counts on this pass
        size_t num_edges = (old_g_iv[old_g_idx + OLD_G_NODE_FROM_COUNT_OFFSET]
                            + old_g_iv[old_g_idx + OLD_G_NODE_TO_COUNT_OFFSET]);
        
        new_g_idx += G_NODE_HEADER_LENGTH + num_edges;
        old_g_idx += G_NODE_HEADER_LENGTH + num_edges * OLD_G_EDGE_LENGTH;
    }
        
    // set up the new rank/selects
    sdsl::util::assign(g_bv_rank, sdsl::rank_support_v<1>(&g_bv));
    sdsl::util::assign(g_bv_select, sdsl::bit_vector::select_1_type(&g_bv));
    
    // now add the edges in using the rank/select operations to translate between offset
    for (size_t old_g_idx = 0, new_g_idx = 0; new_g_idx < g_iv.size(); ) {
        
        size_t num_from = old_g_iv[old_g_idx + OLD_G_NODE_FROM_COUNT_OFFSET];
        size_t num_to = old_g_iv[old_g_idx + OLD_G_NODE_TO_COUNT_OFFSET];
        size_t num_edges = num_to + num_from;
        
#ifdef VERBOSE_DEBUG
        cerr << "updating edges for old index " << old_g_idx << " to new index " << new_g_idx << " with " << num_from << " 'from' edges and " << num_to << " 'to' edges" << endl;
#endif
        
        // recount the edges as either left/right instead of to/from
        size_t num_left = 0, num_right = 0;
        size_t begin = old_g_idx + G_NODE_HEADER_LENGTH;
        size_t end = begin + num_to * OLD_G_EDGE_LENGTH;
        for (size_t i = begin; i < end; i += OLD_G_EDGE_LENGTH) {
            bool from_rev, to_rev;
            orientation_from_old_edge_type(old_g_iv[i + OLD_G_EDGE_TYPE_OFFSET], from_rev, to_rev);
            if (to_rev) {
                ++num_right;
            }
            else {
                ++num_left;
            }
        }
        begin = end;
        end = begin + num_from * OLD_G_EDGE_LENGTH;
        for (size_t i = begin; i < end; i += OLD_G_EDGE_LENGTH) {
            bool from_rev, to_rev;
            orientation_from_old_edge_type(old_g_iv[i + OLD_G_EDGE_TYPE_OFFSET], from_rev, to_rev);
            if (from_rev) {
                ++num_left;
            }
            else {
                ++num_right;
            }
        }
        
#ifdef VERBOSE_DEBUG
        cerr << "translates to " << num_left << " left edges and " << num_right << " right edges" << endl;
#endif
        
        // record the left/right counts
        g_iv[new_g_idx + G_NODE_LEFT_COUNT_OFFSET] = num_left;
        g_iv[new_g_idx + G_NODE_RIGHT_COUNT_OFFSET] = num_right;
        
        // convert the to/from edges into left/right edges
        size_t next_left_idx = new_g_idx + G_NODE_HEADER_LENGTH;
        size_t next_right_idx = next_left_idx + num_left;
        
        // first the to edges
        begin = old_g_idx + G_NODE_HEADER_LENGTH;
        end = begin + num_to * OLD_G_EDGE_LENGTH;
        for (size_t i = begin; i < end; i += OLD_G_EDGE_LENGTH) {
            bool from_rev, to_rev;
            orientation_from_old_edge_type(old_g_iv[i + OLD_G_EDGE_TYPE_OFFSET], from_rev, to_rev);
            
            // apply the relative offset in the old vector
            size_t old_g_nbr_idx = old_g_idx + old_g_iv[i + OLD_G_EDGE_OFFSET_OFFSET];
            // translate that offset into the new vector using rank/select
            size_t new_g_nbr_idx = g_bv_select(old_g_bv_rank(old_g_nbr_idx) + 1);
            
#ifdef VERBOSE_DEBUG
            cerr << "\told 'to' edge with type " << old_g_iv[i + OLD_G_EDGE_TYPE_OFFSET] << " and neighbor at " << old_g_nbr_idx << " now has neighbor " << new_g_nbr_idx << " with to_rev=" << to_rev << " and from_rev=" << from_rev << endl;
#endif
            
            if (!to_rev) {
                g_iv[next_left_idx] = encode_edge(new_g_idx, new_g_nbr_idx, from_rev);
                ++next_left_idx;
            }
            else {
                g_iv[next_right_idx] = encode_edge(new_g_idx, new_g_nbr_idx, !from_rev);
                ++next_right_idx;
            }
        }
        // now the from edges (TODO: repetitive)
        begin = end;
        end = begin + num_from * OLD_G_EDGE_LENGTH;
        for (size_t i = begin; i < end; i += OLD_G_EDGE_LENGTH) {
            bool from_rev, to_rev;
            orientation_from_old_edge_type(old_g_iv[i + OLD_G_EDGE_TYPE_OFFSET], from_rev, to_rev);
            
            // apply the relative offset in the old vector
            size_t old_g_nbr_idx = old_g_idx + old_g_iv[i + OLD_G_EDGE_OFFSET_OFFSET];
            // translate that offset into the new vector using rank/select
            size_t new_g_nbr_idx = g_bv_select(old_g_bv_rank(old_g_nbr_idx) + 1);
            
#ifdef VERBOSE_DEBUG
            cerr << "\told 'from' edge with type " << old_g_iv[i + OLD_G_EDGE_TYPE_OFFSET] << " and neighbor at " << old_g_nbr_idx << " now has neighbor " << new_g_nbr_idx << " with to_rev=" << to_rev << " and from_rev=" << from_rev << endl;
#endif
            
            if (from_rev) {
                g_iv[next_left_idx] = encode_edge(new_g_idx, new_g_nbr_idx, !to_rev);
                ++next_left_idx;
            }
            else {
                g_iv[next_right_idx] = encode_edge(new_g_idx, new_g_nbr_idx, to_rev);
                ++next_right_idx;
            }
        }
        
        new_g_idx += G_NODE_HEADER_LENGTH + num_edges;
        old_g_idx += G_NODE_HEADER_LENGTH + num_edges * OLD_G_EDGE_LENGTH;
    }
        
    sdsl::util::bit_compress(g_iv);
}

size_t XG::edge_index(const edge_t& edge) const {
    // We are going to get our unique indexes by starting with a starting index
    // ofr each node, and then adding to that the number of edges we need to
    // enumerate before we find the one we are looking for, taking all the
    // edges off one side to be enumerated before all the edges off the other.
    // Note that we need to make sure we enumerate edges off the left and right
    // sides of the node in a consistent order, so that they don't end up
    // having colliding indexes.
    
    // Note that resulting indexes will not be anywhere near dense.
    
    edge_t canonical = edge_handle(edge.first, edge.second);
    
    // Get the g index corresponding to the first node's record. We know it
    // owns at least as much g vector space as it has edges.
    // Turns out the g index is just what we packed in the handle's number; no
    // need to do a select here.
    size_t g_start = handlegraph::number_bool_packing::unpack_number(canonical.first);
    
    // Get the 0-based rank so that we know how many node records come before this one.
    size_t prev_nodes = g_bv_rank(g_start);
    
    // Since we know where we are in g_iv, and how many previous nodes there
    // are, and how big a node record is, and how big an edge record is, we
    // can calculate how many edge records came before this node's start. Each
    // edge has two records, and we can't tell how many of those were the
    // canonical records for their edges, but we can at least restrict our
    // non-density in edge index space to a factor of 2 expansion.
    size_t node_start_idx = g_start - (G_NODE_HEADER_LENGTH * prev_nodes);
   
    // Start as the first edge for this node
    size_t idx = node_start_idx;
   
    if (get_is_reverse(canonical.first)) {
        // If we have the first node in locally reverse orientation, add its degree
        // on the locally forward right (or as-presented left) to the index.
        // This avoids collisions between edges attached to different ends of
        // the node.
        idx += get_degree(canonical.first, true);
    }
    
    // Then scan all the edges attached to this end of the node until we find
    // the one we are looking for, and add that count to the index.
    bool not_seen = true;
    
    follow_edges(canonical.first, false, [&](const handle_t& next) {
            // For each edge on the correct side of the node
            
            // Set flag false if we found it
            not_seen = (next != canonical.second);
            
            if (not_seen) {
                // If we didn't find it, increment the index
                ++idx;
            }
            
            return not_seen;
        });
    if (not_seen) {
        // Complain if we don't eventually see the edge.
        throw std::runtime_error("Cound not find index of edge connecting " +
                                 std::to_string(get_id(edge.first)) + " and " + std::to_string(get_id(edge.second)));
    } else {
        // We found it and we have the correct index.
        // TODO: it is 0-based, and vg pack demands 1-based indexes. So we will
        // add 1. See
        // https://github.com/vgteam/libhandlegraph/issues/41#issuecomment-571386849
        return idx + 1;
    }
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

bool XG::follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const {
    
    // unpack the handle
    size_t g = handlegraph::number_bool_packing::unpack_number(handle);
    bool is_reverse = handlegraph::number_bool_packing::unpack_bit(handle);
    
    // get the index range corresponding to the side's edges
    size_t begin, end;
    if (go_left != is_reverse) {
        // left side
        begin = g + G_NODE_HEADER_LENGTH;
        end = begin + g_iv[g + G_NODE_LEFT_COUNT_OFFSET];
    }
    else {
        // right side
        begin = g + G_NODE_HEADER_LENGTH + g_iv[g + G_NODE_LEFT_COUNT_OFFSET];
        end = begin + g_iv[g + G_NODE_RIGHT_COUNT_OFFSET];
    }
    
    // construct the neighboring handles and execute the iteratee on them
    bool keep_going = true;
    for (size_t i = begin; i < end && keep_going; ++i) {
        uint64_t raw_edge = g_iv[i];
        handle_t next = number_bool_packing::pack(g + edge_relative_offset(raw_edge),
                                                  edge_orientation(raw_edge) != is_reverse);
        keep_going = iteratee(next);
    }
    return keep_going;
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
                    
                    // This record is the header plus all the edge records it contains.
                    // Decode the entry size in the same thread doing the iteration.
                    g += G_NODE_HEADER_LENGTH + g_iv[g + G_NODE_LEFT_COUNT_OFFSET] + g_iv[g + G_NODE_RIGHT_COUNT_OFFSET];
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
            
            // This record is the header plus all the edge records it contains.
            // Decode the entry size in the same thread doing the iteration.
            g += G_NODE_HEADER_LENGTH + g_iv[g + G_NODE_LEFT_COUNT_OFFSET] + g_iv[g + G_NODE_RIGHT_COUNT_OFFSET];
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
    std::string query = path_name_csa_delim + path_name + path_name_csa_delim;
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
    size_t start = pn_bv_select(rank)+1; // step past '$'
    size_t end = pn_bv_select(rank+1);
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

size_t XG::get_edge_count() const {
    return this->edge_count;
}

size_t XG::get_total_length() const {
    return this->seq_length;
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
    // We expect input positions as 1-based for now.
    // See https://github.com/vgteam/libhandlegraph/issues/41
   
    // We have a rank primitive that gets us the 1s *before* a position.
    // We want all the bases *at or after* a 1 to give us the same value.
    // So we get the 1s before (offset - 1) + 1, which is the same as at or
    // after offset - 1, which is the 0-based version of offset.
    return rank_to_id(s_bv_rank(offset));
}

size_t XG::node_vector_offset(const nid_t& id) const {
    // We produce offsets as 0-based for now.
    // See https://github.com/vgteam/libhandlegraph/issues/41
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

int get_thread_count(void) {
    int thread_count = 1;
#pragma omp parallel
    {
#pragma omp master
        thread_count = omp_get_num_threads();
    }
    return thread_count;
}

}
