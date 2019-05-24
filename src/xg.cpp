#include "xg.hpp"

#include "mmmultimap.hpp"

#include <bitset>
#include <arpa/inet.h>

#include <handlegraph/util.hpp>

#include "gfakluge.hpp"

#define VERBOSE_DEBUG
//#define debug_algorithms
//#define debug_component_index

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
    
    if (file_version > MAX_INPUT_VERSION) {
        // This XG file is from the distant future.
        throw XGFormatError("XG index file version " + std::to_string(file_version) +
                            " is too new to interpret (max version = " + std::to_string(MAX_INPUT_VERSION) + ")");
    }
    
    try {
        
        ////////////////////////////////////////////////////////////////////////
        // DO NOT CHANGE THIS CODE without creating a new XG version:
        // 1. Increment OUTPUT_VERSION to a new integer.
        // 2. Change the serialization code.
        // 3. Add a new case here (or enhance an existing case) with new deserialization code.
        // 4. Up MAX_INPUT_VERSION to allow your new version to be read.
        ////////////////////////////////////////////////////////////////////////
        switch (file_version) {
        
        case 5: // Fall through
        case 6:
        case 7:
        case 8:
        case 9:
            std::cerr << "warning:[XG] Loading an out-of-date XG format."
                      << "For better performance over repeated loads, consider recreating this XG index." << std::endl;
            // Fall through
        case 10:
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
                    path->load(in);
                    paths.push_back(path);
                }
                np_iv.load(in);
                np_bv.load(in);
                np_bv_rank.load(in, &np_bv);
                np_bv_select.load(in, &np_bv);
                
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
    sdsl::read_member(min_node_id, in);
    ids.load(in);
    directions.load(in);
    ranks.load(in);
    positions.load(in);
    offsets.load(in);
    offsets_rank.load(in, &offsets);
    offsets_select.load(in, &offsets);
    sdsl::read_member(is_circular, in);    
}

size_t XGPath::serialize(std::ostream& out,
                         sdsl::structure_tree_node* v,
                         std::string name) const {
    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_t written = 0;
    written += sdsl::write_member(min_node_id, out, child, "min_node_id" + name);
    written += ids.serialize(out, child, "path_node_ids_" + name);
    written += directions.serialize(out, child, "path_node_directions_" + name);
    written += ranks.serialize(out, child, "path_mapping_ranks_" + name);
    written += positions.serialize(out, child, "path_node_offsets_" + name);
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
               XG& graph,
               size_t* unique_member_count_out) {

    // The circularity flag is just a normal bool
    this->is_circular = is_circular;

    // node ids, the literal path
    sdsl::int_vector<> ids_iv;
    sdsl::util::assign(ids_iv, sdsl::int_vector<>(path.size()));
    // directions of traversal (typically forward, but we allow backwards
    sdsl::bit_vector directions_bv;
    sdsl::util::assign(directions_bv, sdsl::bit_vector(path.size()));
    // node positions in path
    sdsl::util::assign(positions, sdsl::int_vector<>(path.size()));
    // mapping ranks in path
    sdsl::util::assign(ranks, sdsl::int_vector<>(path.size()));

    size_t path_off = 0;
    size_t members_off = 0;
    size_t positions_off = 0;
    size_t path_length = 0;
    min_node_id = std::numeric_limits<int64_t>::max();

    // determine min node id
    for (size_t i = 0; i < path.size(); ++i) {
        auto node_id = graph.get_id(path[i]);
        min_node_id = std::min(min_node_id, node_id);
    }
    // determine total length and record node ids
    for (size_t i = 0; i < path.size(); ++i) {
        const handle_t& handle = path[i];
        nid_t node_id = graph.get_id(handle);
        path_length += graph.get_length(handle);
        ids_iv[i] = local_id(node_id);
        // we will explode if the node isn't in the graph
    }

    // make the bitvector for path offsets
    sdsl::util::assign(offsets, sdsl::bit_vector(path_length));
    std::set<int64_t> uniq_nodes;
    //cerr << "path " << path_name << " has " << path.size() << endl;
    for (size_t i = 0; i < path.size(); ++i) {
        //cerr << i << endl;
        auto& handle = path[i];
        auto node_id = graph.get_id(handle);
        bool is_reverse = graph.get_is_reverse(handle);
        // record direction of passage through node
        directions_bv[i] = is_reverse;
        // and the external rank of the mapping
        ranks[i] = number_bool_packing::unpack_number(handle);
        // we've seen another entity
        uniq_nodes.insert(node_id);
        // and record node offset in path
        positions[positions_off++] = path_off;
        // record position of node
        offsets[path_off] = 1;
        // and update the offset counter
        path_off += graph.get_length(handle);
    }
    //cerr << uniq_nodes.size() << " vs " << path.size() << endl;
    if(unique_member_count_out) {
        // set member count as the unique entities that are in the path
        // We don't need it but our caller might
        *unique_member_count_out = uniq_nodes.size();
    }
    // and traversal information
    sdsl::util::assign(directions, sdsl::sd_vector<>(directions_bv));
    // handle entity lookup structure (wavelet tree)
    sdsl::util::bit_compress(ids_iv);
    sdsl::construct_im(ids, ids_iv);
    // bit compress the positional offset info
    sdsl::util::bit_compress(positions);
    // bit compress mapping ranks
    sdsl::util::bit_compress(ranks);

    // and set up rank/select dictionary on them
    sdsl::util::assign(offsets_rank, sdsl::rank_support_v<1>(&offsets));
    sdsl::util::assign(offsets_select, sdsl::bit_vector::select_1_type(&offsets));
}

nid_t XGPath::node(size_t offset) const {
    return external_id(ids[offset]);
}
    
nid_t XGPath::node_at_position(size_t pos) const {
    return node(offset_at_position(pos));
}

size_t XGPath::offset_at_position(size_t pos) const {
    return offsets_rank(pos+1)-1;
}

bool XGPath::is_reverse(size_t offset) const {
    return directions[offset];
}

nid_t XGPath::local_id(nid_t id) const {
    if (id < min_node_id) {
        return std::numeric_limits<int64_t>::max();
    } else {
        return id-min_node_id+1;
    }
}

nid_t XGPath::external_id(nid_t id) const {
    return id+min_node_id-1;
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
    int32_t version_buffer = htonl(OUTPUT_VERSION);
    out.write((char*) &version_buffer, sizeof(version_buffer));
    written += sizeof(version_buffer) / sizeof(char);
    
    ////////////////////////////////////////////////////////////////////////
    // DO NOT CHANGE THIS CODE without creating a new XG version:
    // 1. Increment OUTPUT_VERSION to a new integer.
    // 2. Add your new serialization code.
    // 3. Add a new case for your new version to XG::load()
    // 4. Up MAX_INPUT_VERSION to allow your new version to be read.
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
        paths_written += path->serialize(out, paths_child, "path:" + get_path_name(as_path_handle(i)));
    }
    
    paths_written += np_iv.serialize(out, paths_child, "node_path_mapping");
    paths_written += np_bv.serialize(out, paths_child, "node_path_mapping_starts");
    paths_written += np_bv_rank.serialize(out, paths_child, "node_path_mapping_starts_rank");
    paths_written += np_bv_select.serialize(out, paths_child, "node_path_mapping_starts_select");
    
    sdsl::structure_tree::add_size(paths_child, paths_written);
    written += paths_written;

    sdsl::structure_tree::add_size(child, written);
    return written;
    
}

void XG::from_gfa(const std::string& gfa_filename, std::string basename) {
    if (basename.empty()) {
        basename = gfa_filename;
    }
    char* filename = (char*)gfa_filename.c_str();
    gfak::GFAKluge gfa;
    node_count = 0;
    seq_length = 0;
    edge_count = 0;
    path_count = 0;
    // get information about graph size and id ranges
    gfa.for_each_sequence_line_in_file(filename, [&](gfak::sequence_elem s) {
            nid_t id = std::stol(s.name);
            // min id starts at 0
            if (min_id) {
                min_id = std::min(min_id, id);
            } else {
                min_id = id;
            }
            max_id = std::max(max_id, id);
            seq_length += s.sequence.size();
            ++node_count;
        });
    // edge count
    gfa.for_each_edge_line_in_file(filename, [&](gfak::edge_elem e) {
            if (e.source_name.empty()) return;
            ++edge_count;
        });
    // path count
    gfa.for_each_path_element_in_file(filename, [&](const std::string& path_name_raw, const std::string& node_id, bool is_rev, const std::string& cigar) {
            ++path_count;
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
    sdsl::util::assign(s_bv, sdsl::bit_vector(seq_length));
    sdsl::util::assign(i_iv, sdsl::int_vector<>(node_count));
    sdsl::util::assign(r_iv, sdsl::int_vector<>(max_id-min_id+1)); // note: possibly discontiguous
    
    // for each node in the sequence
    // concatenate the labels into the s_iv
#ifdef VERBOSE_DEBUG
    cerr << "storing node labels" << endl;
#endif
    size_t r = 1;    
    // first make i_iv and r_iv
    gfa.for_each_sequence_line_in_file(filename, [&](gfak::sequence_elem s) {
            nid_t id = std::stol(s.name);
            i_iv[r-1] = id;
            // store ids to rank mapping
            r_iv[id-min_id] = r;
            ++r;
        });
    sdsl::util::bit_compress(i_iv);
    sdsl::util::bit_compress(r_iv);

    // then make s_bv and s_iv
    gfa.for_each_sequence_line_in_file(filename, [&](gfak::sequence_elem s) {
            nid_t id = std::stol(s.name);
            size_t i = r_iv[id-min_id]-1;
            s_bv[i] = 1; // record node start
            for (auto c : s.sequence) {
                s_iv[i++] = dna3bit(c); // store sequence
            }
        });

    // to label the paths we'll need to compress and index our vectors
    sdsl::util::bit_compress(s_iv);
    sdsl::util::assign(s_bv_rank, sdsl::rank_support_v<1>(&s_bv));
    sdsl::util::assign(s_bv_select, sdsl::bit_vector::select_1_type(&s_bv));
    
    // now that we've set up our sequence indexes, we can build the locally traversable graph storage

    // first, we need to collect the edges for each node
    // we use the mmmultimap here to reduce in-memory costs to a minimum
    std::string edge_f_t_idx = basename + ".from_to.mm";
    std::string edge_t_f_idx = basename + ".to_from.mm";
    multimap<uint64_t, uint64_t> edge_from_to_mm(edge_f_t_idx);
    multimap<uint64_t, uint64_t> edge_to_from_mm(edge_t_f_idx);
    gfa.for_each_edge_line_in_file(filename, [&](gfak::edge_elem e) {
            if (e.source_name.empty()) return;
            uint64_t from_rank = r_iv[std::stol(e.source_name)-min_id];
            handle_t from_handle = number_bool_packing::pack(from_rank, !e.source_orientation_forward);
            uint64_t to_rank = r_iv[std::stol(e.source_name)-min_id];
            handle_t to_handle = number_bool_packing::pack(to_rank, !e.sink_orientation_forward);
            edge_from_to_mm.append(as_integer(from_handle), as_integer(to_handle));
            edge_to_from_mm.append(as_integer(to_handle), as_integer(from_handle));
        });
    handle_t max_handle = number_bool_packing::pack(max_id, true);
    edge_from_to_mm.index(as_integer(max_handle));
    edge_to_from_mm.index(as_integer(max_handle));

    // calculate g_iv size
    size_t g_iv_size =
        node_count * G_NODE_HEADER_LENGTH // record headers
        + edge_count * 2 * G_EDGE_LENGTH; // edges (stored twice)
    sdsl::util::assign(g_iv, sdsl::int_vector<>(g_iv_size));
    sdsl::util::assign(g_bv, sdsl::bit_vector(g_iv_size));

    auto temp_get_handle = [&](const nid_t& id, bool orientation) {
        uint64_t handle_rank = r_iv[id-min_id];
        return number_bool_packing::pack(handle_rank, orientation);
    };
    auto temp_get_id = [&](const handle_t& h) {
        return i_iv[number_bool_packing::unpack_number(h)];
    };
    
    int64_t g = 0; // pointer into g_iv and g_bv
    for (int64_t i = 0; i < node_count; ++i) {
        nid_t id = i_iv[i];
        //std::cerr << "id is " << id << std::endl;
        handle_t handle = temp_get_handle(id, false);
        g_bv[g] = 1; // mark record start for later query
        g_iv[g++] = id;
        g_iv[g++] = node_start(id);
        g_iv[g++] = get_sequence(handle).size();
        size_t to_edge_count = 0;
        size_t from_edge_count = 0;
        size_t to_edge_count_idx = g++;
        size_t from_edge_count_idx = g++;
        // write the edges in id-based format
        // we will next convert these into relative format
        for (auto orientation : { false, true }) {
            handle_t to = temp_get_handle(id, orientation);
            edge_to_from_mm.for_unique_values_of(as_integer(to), [&](const uint64_t& _from) {
                    handle_t from = as_handle(_from);
                    g_iv[g++] = temp_get_id(from);
                    g_iv[g++] = edge_type(from, to);
                    ++to_edge_count;
                });
        }
        g_iv[to_edge_count_idx] = to_edge_count;
        for (auto orientation : { false, true }) {
            handle_t from = temp_get_handle(id, orientation);
            edge_to_from_mm.for_unique_values_of(as_integer(from), [&](const uint64_t& _to) {
                    handle_t to = as_handle(_to);
                    g_iv[g++] = temp_get_id(to);
                    g_iv[g++] = edge_type(from, to);
                    ++from_edge_count;
                });
        }
        g_iv[from_edge_count_idx] = from_edge_count;
    }

    // cleanup our mmmultimap
    std::remove(edge_f_t_idx.c_str());
    std::remove(edge_t_f_idx.c_str());
    
    // set up rank and select supports on g_bv so we can locate nodes in g_iv
    sdsl::util::assign(g_bv_rank, sdsl::rank_support_v<1>(&g_bv));
    sdsl::util::assign(g_bv_select, sdsl::bit_vector::select_1_type(&g_bv));

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
    string path_names;
    size_t path_node_count = 0; // count of node path memberships

    std::string curr_path_name;
    std::vector<handle_t> curr_path_steps;
    size_t curr_node_count = 0;
    bool curr_is_circular = false; // TODO, use TP:Z:circular tag... we'll have to fish this out of the file

    auto build_accumulated_path = [&](void) {
        size_t unique_member_count = 0;
        path_names += start_marker + curr_path_name + end_marker;
        XGPath* path = new XGPath(curr_path_name, curr_path_steps,
                                  curr_is_circular,
                                  *this, &unique_member_count);
        paths.push_back(path);
        path_node_count += unique_member_count;
    };

    // todo ... is it circular?
    // might make sense to scan the file for this
    gfa.for_each_path_element_in_file(filename, [&](const std::string& path_name_raw, const std::string& node_id, bool is_rev, const std::string& cigar) {
            std::string path_name = path_name_raw;
            path_name.erase(std::remove_if(path_name.begin(), path_name.end(), [](char c) { return std::isspace(c); }), path_name.end());
            if (path_name != curr_path_name && !curr_path_name.empty()) {
                // build the last path we've accumulated
                build_accumulated_path();
                curr_path_steps.clear();
                curr_path_name = path_name;
            }
            curr_path_steps.push_back(get_handle(stol(node_id), is_rev));
        });
    // build the last path
    build_accumulated_path();
    curr_path_steps.clear();

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
    string path_name_file = "@pathnames.iv";
    sdsl::store_to_file((const char*)path_names.c_str(), path_name_file);
    sdsl::construct(pn_csa, path_name_file, 1);

    // node -> paths
    sdsl::util::assign(np_iv, sdsl::int_vector<>(path_node_count+node_count));
    sdsl::util::assign(np_bv, sdsl::bit_vector(path_node_count+node_count));
    size_t np_off = 0;
    for (size_t i = 0; i < node_count; ++i) {
        np_bv[np_off] = 1;
        np_iv[np_off] = 0; // null so we can detect entities with no path membership
        ++np_off;
        nid_t id = rank_to_id(i+1);
        for (size_t j = 1; j <= paths.size(); ++j) {
            if (node_occs_in_path(id, as_path_handle(j)) > 0) {
                np_iv[np_off++] = j;
            }
        }
    }

    sdsl::util::bit_compress(np_iv);
    //cerr << ep_off << " " << path_entities << " " << entity_count << endl;
    assert(np_off <= path_node_count+node_count);
    sdsl::util::assign(np_bv_rank, sdsl::rank_support_v<1>(&np_bv));
    sdsl::util::assign(np_bv_select, sdsl::bit_vector::select_1_type(&np_bv));

    // memoize which paths co-occur on connected components
    // todo... implement this again??
    //index_component_path_sets();

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
            for (size_t j = 0; j + 1 < path->ids.size(); j++) {
                cerr << path->node(j) << " ";
            }
            if (path->ids.size() > 0) {
                cerr << path->node(path->ids.size() - 1);
            }
            cerr << endl;
            cerr << path->ranks << endl;
            cerr << path->directions << endl;
            cerr << path->positions << endl;
            cerr << path->offsets << endl;
        }
        cerr << np_bv << endl;
        cerr << np_iv << endl;
        
    }

    // TODO check if we've made a valid graph
    //
}
    
char XG::pos_char(nid_t id, bool is_rev, size_t off) const {
    assert(off < get_length(id));
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

size_t XG::id_to_rank(const nid_t& id) const {
    size_t x = id-min_id;
    if (x < 0 || x >= r_iv.size()) return 0;
    return r_iv[id-min_id];
}

nid_t XG::rank_to_id(const size_t& rank) const {
    if(rank == 0) {
        cerr << "[xg] error: Request for id of rank 0" << endl;
        assert(false);
    }
    if(rank > node_count) {
        cerr << "[xg] error: Request for id of rank " << rank << "/" << node_count << endl;
        assert(false);
    }
    return g_iv[g_bv_select(rank)];
}

int XG::edge_type(const handle_t& from, const handle_t& to) const {
    if (get_is_reverse(from) && get_is_reverse(to)) {
        return 4;
    } else if (get_is_reverse(from)) {
        return 3;
    } else if (get_is_reverse(to)) {
        return 2;
    } else {
        return 1;
    }
}

handle_t XG::get_handle(const nid_t& node_id, bool is_reverse) const {
    // Handles will be g vector index with is_reverse in the low bit
    
    // Where in the g vector do we need to be
    uint64_t g = g_bv_select(id_to_rank(node_id));
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

bool XG::edge_filter(int type, bool is_to, bool want_left, bool is_reverse) const {
    // Return true if we want an edge of the given type, where we are the from
    // or to node (according to is_to), when we are looking off the right or
    // left side of the node (according to want_left), and when the node is
    // forward or reverse (accoridng to is_reverse).
    
    // Edge type encoding:
    // 1: end to start
    // 2: end to end
    // 3: start to start
    // 4: start to end
    
    // First compute what we want looking off the right of a node in the forward direction.
    bool wanted = !is_to && (type == 1 || type == 2) || is_to && (type == 2 || type == 4);
    
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
        assert(type >= 0);
        assert(type <= 3);
        
        if (edge_filter(type, is_to, want_left, is_reverse)) {
            
            // What's the offset to the other node?
            int64_t offset = g_iv[start + i * G_EDGE_LENGTH + G_EDGE_OFFSET_OFFSET];
            
            // Make sure we haven't gone off the rails into non-edge data.
            assert((int64_t) g + offset >= 0);
            assert(g + offset < g_iv.size());
            
            // Should we invert?
            // We only invert if we cross an end to end edge. Or a start to start edge
            bool new_reverse = is_reverse != (type == 2 || type == 3);
            
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
            bool new_reverse = is_reverse != (type == 2 || type == 3);
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
                        // Run the iteratee
                        if (!iteratee(handle)) {
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
            
            // Run the iteratee in-line
            if (!iteratee(handle)) {
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
        cerr << "multiple hits for " << query << endl;
        assert(false);
    }
    if(occs.size() == 0) {
        // This path does not exist. Give back 0, which can never be a real path
        // rank.
        return as_path_handle(0);
    }
    return as_path_handle(pn_bv_rank(occs[0])+1); // step past '#'
}
    
std::string XG::get_path_name(const path_handle_t& path_handle) const {
    return path_name(as_integer(path_handle));
}

bool XG::get_is_circular(const path_handle_t& path_handle) const {
    return paths[as_integer(path_handle) - 1]->is_circular;
}

size_t XG::get_step_count(const path_handle_t& path_handle) const {
    return paths[as_integer(path_handle) - 1]->ids.size();
}
    
size_t XG::get_path_length(const path_handle_t& path_handle) const {
    return paths[as_integer(path_handle) - 1]->offsets.size();
}

handle_t XG::get_handle_of_step(const step_handle_t& step_handle) const {
    const auto& xgpath = *paths[as_integer(get_path_handle_of_step(step_handle)) - 1];
    size_t idx = as_integers(step_handle)[1];
    return get_handle(xgpath.node(idx), xgpath.is_reverse(idx));
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
    
    std::vector<std::pair<path_handle_t, std::vector<std::pair<size_t, bool>>>> oriented_paths = oriented_paths_of_node(get_id(handle));
    
    for (const std::pair<path_handle_t, std::vector<std::pair<size_t, bool>>>& path_occs : oriented_paths) {
        for (const std::pair<size_t, bool>& oriented_occ : path_occs.second) {
            step_handle_t step_handle;
            as_integers(step_handle)[0] = as_integer(path_occs.first);
            as_integers(step_handle)[1] = oriented_occ.first;
            if (!iteratee(step_handle)) {
                return false;
            }
        }
    }
    
    return true;
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

int64_t XG::node_at_seq_pos(const size_t& pos) const {
    return rank_to_id(s_bv_rank(pos));
}

size_t XG::node_start(const nid_t& id) const {
    return s_bv_select(id_to_rank(id));
}

size_t XG::node_graph_idx(const nid_t& id) const {
    return g_bv_select(id_to_rank(id));
}

const XGPath& XG::get_path(const std::string& name) const {
    return *paths[as_integer(get_path_handle(name))-1];
}

std::vector<size_t> XG::path_ranks_by_prefix(const std::string& prefix) const {
    // find the name in the csa
    std::string query = start_marker + prefix;
    auto occs = locate(pn_csa, query);
    std::vector<size_t> ranks;
    for (size_t i = 0; i < occs.size(); ++i) {
        ranks.push_back(pn_bv_rank(occs[i])+1); // step past '#'
    }
    return ranks;
}

std::vector<std::string> XG::path_names_by_prefix(const std::string& prefix) const {
    std::vector<std::string> names;
    for (auto& rank : path_ranks_by_prefix(prefix)) {
        names.push_back(path_name(rank));
    }
    return names;
}

std::string XG::path_name(const size_t& rank) const {
    size_t start = pn_bv_select(rank)+1; // step past '#'
    size_t end = rank == path_count ? pn_iv.size() : pn_bv_select(rank+1);
    end -= 1;  // step before '$'
    string name; name.resize(end-start);
    for (size_t i = start; i < end; ++i) {
        name[i-start] = pn_iv[i];
    }
    return name;
}

bool XG::path_contains_node(const std::string& name, const nid_t& id) const {
    return node_occs_in_path(id, name) > 0;
}

std::vector<path_handle_t> XG::paths_of_node(const nid_t& id) const {
    auto rank = id_to_rank(id);
    if (rank == 0) {
        throw runtime_error("Tried to get paths of nonexistent node " + to_string(id));
    }
    size_t off = np_bv_select(rank);
    assert(np_bv[off++]);
    std::vector<path_handle_t> path_handles;
    while (off < np_bv.size() ? np_bv[off] == 0 : false) {
        path_handles.push_back(as_path_handle(np_iv[off++]));
    }
    return path_handles;
}

std::vector<std::pair<path_handle_t, std::vector<std::pair<size_t, bool>>>> XG::oriented_paths_of_node(const nid_t& id) const {
    std::vector<path_handle_t> node_paths = paths_of_node(id);
    return oriented_occurrences_on_paths(id, node_paths);
}
    
std::vector<std::pair<size_t, bool>> XG::oriented_occurrences_on_path(const nid_t& id, const path_handle_t& path) const {
    std::vector<std::pair<size_t, bool>> occurrences;
    for (size_t i : node_ranks_in_path(id, path)) {
        occurrences.emplace_back(i, paths[as_integer(path)-1]->directions[i]);
    }
    return occurrences;
}
    
std::vector<std::pair<path_handle_t, std::vector<std::pair<size_t, bool>>>>
XG::oriented_occurrences_on_paths(const nid_t& id, std::vector<path_handle_t>& paths) const {
    std::vector<std::pair<path_handle_t, std::vector<std::pair<size_t, bool>>>> path_occurrences;
    for (auto& path_handle : paths) {
        path_occurrences.emplace_back(path_handle, oriented_occurrences_on_path(id, path_handle));
    }
    return path_occurrences;
}

std::pair<pos_t, int64_t> XG::next_path_position(const pos_t& pos, const int64_t& max_search) const {
    handle_t h_fwd = get_handle(id(pos), is_rev(pos));
    handle_t h_rev = get_handle(id(pos), !is_rev(pos));
    int64_t fwd_seen = offset(pos);
    int64_t rev_seen = get_length(h_fwd) - offset(pos);
    std::pair<pos_t, int64_t> fwd_next = make_pair(make_pos_t(0,false,0), numeric_limits<int64_t>::max());
    std::pair<pos_t, int64_t> rev_next = make_pair(make_pos_t(0,false,0), numeric_limits<int64_t>::max());
    follow_edges(h_fwd, false, [&](const handle_t& n) {
            nid_t id = get_id(n);
            if (!paths_of_node(id).empty()) {
                fwd_next = make_pair(make_pos_t(id, get_is_reverse(n), 0), fwd_seen);
                return false;
            } else {
                fwd_seen += get_length(n);
                return fwd_seen < max_search;
            }
        });
    follow_edges(h_rev, false, [&](const handle_t& n) {
            nid_t id = get_id(n);
            if (!paths_of_node(id).empty()) {
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
            std::vector<path_handle_t> path_ids = paths_of_node(id);
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
            auto name = get_path_name(i.first);
            for (auto& p1 : position_in_path(near1.first, name)) {
                for (auto& p2 : position_in_path(near2.first, name)) {
                    int64_t distance = abs((int64_t)p1 - (int64_t)p2);
                    min_distance = min(distance, min_distance);
                }
            }
        }
    }
    return min_distance;
}

void XG::for_path_range(const std::string& name, int64_t start, int64_t stop,
                        std::function<void(int64_t, bool)> lambda, bool is_rev) const {
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
        lambda(path.node(i), path.is_reverse(i));
    }
}

size_t XG::node_occs_in_path(const nid_t& id, const std::string& name) const {
    return node_occs_in_path(id, get_path_handle(name));
}

size_t XG::node_occs_in_path(const nid_t& id, const path_handle_t& path) const {
    size_t p = as_integer(path)-1;
    auto& pi_wt = paths[p]->ids;
    int64_t local_id = paths[p]->local_id(id);
    return pi_wt.rank(pi_wt.size(), local_id);
}

std::vector<size_t> XG::node_ranks_in_path(const nid_t& id, const std::string& name) const {
    return node_ranks_in_path(id, get_path_handle(name));
}

std::vector<size_t> XG::node_ranks_in_path(const nid_t& id, const path_handle_t& path) const {
    std::vector<size_t> ranks;
    size_t p = as_integer(path)-1;
    size_t occs = node_occs_in_path(id, path);
    int64_t local_id = paths[p]->local_id(id);
    for (size_t i = 1; i <= occs; ++i) {
        ranks.push_back(paths[p]->ids.select(i, local_id));
    }
    return ranks;
}

std::vector<size_t> XG::position_in_path(const nid_t& id, const std::string& name) const {
    return position_in_path(id, get_path_handle(name));
}

std::vector<size_t> XG::position_in_path(const nid_t& id, const path_handle_t& path) const {
    auto& _path = *paths[as_integer(path)-1];
    std::vector<size_t> pos_in_path;
    for (auto i : node_ranks_in_path(id, path)) {
        pos_in_path.push_back(_path.positions[i]);
    }
    return pos_in_path;
}

std::map<std::string, std::vector<size_t> > XG::position_in_paths(const nid_t& id, bool is_rev, const size_t& offset) const {
    std::map<std::string, std::vector<size_t> > positions;
    for (auto& prank : paths_of_node(id)) {
        auto& path = *paths[as_integer(prank)-1];
        auto& pos_in_path = positions[get_path_name(prank)];
        for (auto i : node_ranks_in_path(id, prank)) {
            size_t pos = offset + (is_rev ?
                                   // Account for the reverse-strand offset
                                   get_path_length(prank) - path.positions[i] - get_length(get_handle(id, false))
                                   : path.positions[i]);
            pos_in_path.push_back(pos);
        }
    }
    return positions;
}

std::map<std::string, std::vector<std::pair<size_t, bool> > > XG::offsets_in_paths(pos_t pos) const {
    std::map<std::string, std::vector<std::pair<size_t, bool> > > positions;
    nid_t node_id = id(pos);
    for (auto& prank : paths_of_node(node_id)) {
        auto& path = *paths[as_integer(prank)-1];
        auto& pos_in_path = positions[get_path_name(prank)];
        for (auto i : node_ranks_in_path(node_id, prank)) {
            // relative direction to this traversal
            bool dir = path.directions[i] != is_rev(pos);
            // Make sure to interpret the pos_t offset on the correct strand.
            // Normalize to a forward strand offset.
            size_t node_forward_strand_offset = is_rev(pos) ? (get_length(get_handle(node_id, false)) - offset(pos) - 1) : offset(pos);
            // Then go forward or backward along the path as appropriate. If
            // the node is on the path in reverse we have where its end landed
            // and have to flip the forward strand offset around again.
            size_t off = path.positions[i] + (path.directions[i] ?
                                              (get_length(get_handle(node_id, false)) - node_forward_strand_offset - 1) :
                                              node_forward_strand_offset);
            pos_in_path.push_back(make_pair(off, dir));
        }
    }
    return positions;
}

std::map<std::string, std::vector<std::pair<size_t, bool> > > XG::nearest_offsets_in_paths(const pos_t& pos, int64_t max_search) const {
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
        std::map<std::string, std::vector<std::pair<size_t, bool> > > empty;
        return empty;
    }
}

std::map<std::string, std::vector<size_t> > XG::distance_in_paths(const nid_t& id1, bool is_rev1, const size_t& offset1,
                                                                  const nid_t& id2, bool is_rev2, const size_t& offset2) const {
    auto pos1 = position_in_paths(id1, is_rev1, offset1);
    auto pos2 = position_in_paths(id2, is_rev2, offset2);
    std::map<std::string, std::vector<size_t> > dist;
    // position in a path is undefined in inversion
    if (is_rev1 != is_rev2) {
        return dist;
    }
    for (auto& c1 : pos1) {
        auto c2 = pos2.find(c1.first);
        if (c2 != pos2.end()) {
            auto& d = dist[c1.first];
            // distances are a cross of the points
            for (auto o1 : c1.second) {
                for (auto o2 : c2->second) {
                    d.push_back(o1-o2);
                }
            }
        }
    }
    return dist;
}

int64_t XG::min_distance_in_paths(const nid_t& id1, bool is_rev1, const size_t& offset1,
                                  const nid_t& id2, bool is_rev2, const size_t& offset2) const {
    auto dist = distance_in_paths(id1, is_rev1, offset1,
                                  id2, is_rev2, offset2);
    size_t min_dist = std::numeric_limits<size_t>::max();
    for (auto& c : dist) {
        for (auto& o : c.second) {
            if (o <  min_dist) {
                min_dist = o;
            }
        }
    }
    return min_dist;
}

nid_t XG::node_at_path_position(const std::string& name, size_t pos) const {
    return paths[as_integer(get_path_handle(name))]->node_at_position(pos);
}

size_t XG::node_start_at_path_position(const std::string& name, size_t pos) const {
    size_t p = as_integer(get_path_handle(name))-1;
    size_t position_rank = paths[p]->offsets_rank(pos+1);
    return paths[p]->offsets_select(position_rank);
}

pos_t XG::graph_pos_at_path_position(const std::string& name, size_t path_pos) const {
    auto& path = get_path(name);
    path_pos = min((size_t)path.offsets.size()-1, path_pos);
    size_t trav_idx = path.offsets_rank(path_pos+1)-1;
    // Get the offset along the node in its path direction.
    // If the node is forward along the path, we get the forward strand offset on the node, and return a forward pos_t.
    // If the node is backward along the path, we get the reverse strand offset automatically, and return a reverse pos_t.
    int64_t offset = path_pos - path.positions[trav_idx];
    nid_t node_id = path.node(trav_idx);
    bool is_rev = path.directions[trav_idx];
    return make_pos_t(node_id, is_rev, offset);
}

}
