#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

struct maf_row_t {
    std::string path_name;
    uint64_t record_start = 0;
    uint64_t seq_size = 0;
    bool is_reversed = false;
    uint64_t path_length = 0;
    std::string aligned_seq;
};

struct maf_partial_row_t {
    uint64_t record_start = 0;
    uint64_t seq_size = 0;
    bool is_reversed = false;
    uint64_t path_length = 0;
    std::string aligned_seq;
};

struct maf_t {
    std::vector<uint64_t> block_ids;
    std::map<std::string, std::vector<maf_partial_row_t>> rows;
    std::vector<std::pair<std::string, maf_partial_row_t>> consensus_rows;
};

template<typename T>
void clear_vector(std::vector<T>& vec) {
    vec.clear();
    vec.shrink_to_fit();
    std::vector<T>().swap(vec);
}
void clear_string(std::string str){
    str.clear();
    str.shrink_to_fit();
    std::string().swap(str);
}

void write_maf_rows(std::ofstream &out, const std::map<std::string, std::vector<maf_partial_row_t>>& maf) {
    // determine output widths for everything
    size_t max_src_length = 0;
    size_t max_start_length = 0;
    size_t max_seq_size_length = 0;
    size_t max_is_rev_length = 1;
    size_t max_src_size_length = 0;
    size_t max_text_length = 0;
    for (auto& path_to_maf_rows : maf) {
        for (auto &row : path_to_maf_rows.second) {
            max_src_length = std::max(max_src_length, path_to_maf_rows.first.size());
            max_start_length = std::max(max_start_length, std::to_string(row.record_start).size());
            max_seq_size_length = std::max(max_seq_size_length, std::to_string(row.seq_size).size());
            max_src_size_length = std::max(max_src_size_length, std::to_string(row.path_length).size());
            max_text_length = std::max(max_text_length, row.aligned_seq.size());
        }
    }

    for (auto& path_to_maf_rows : maf) {
        // write and pad them
        for (auto &row : path_to_maf_rows.second) {
            out << "s "
                << path_to_maf_rows.first << std::string(max_src_length - path_to_maf_rows.first.size(), ' ')
                << std::setw(max_start_length+1) << row.record_start
                << std::setw(max_seq_size_length+1) << row.seq_size
                << std::setw(max_is_rev_length+1) << (row.is_reversed ? "-" : "+")
                << std::setw(max_src_size_length+1) << row.path_length
                << " " << row.aligned_seq
                << "\n";
        }
    }
    out << std::endl;
}
