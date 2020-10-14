#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

struct maf_row_t {
    std::string path_name;
    uint64_t record_start;
    uint64_t seq_size;
    bool is_reversed;
    uint64_t path_length;
    std::string aligned_seq;
};

struct maf_partial_row_t {
    uint64_t record_start;
    uint64_t seq_size;
    bool is_reversed;
    uint64_t path_length;
    std::string aligned_seq;
};

struct maf_t {
    std::vector<uint64_t> field_blocks;
    std::map<std::string, maf_partial_row_t> rows;
    std::vector<std::pair<std::string, maf_partial_row_t>> consensus_rows;
};

void write_maf_rows(std::ofstream &out, const std::vector<maf_row_t>& rows) {
    // determine output widths for everything
    size_t max_src_length = 0;
    size_t max_start_length = 0;
    size_t max_seq_size_length = 0;
    size_t max_is_rev_length = 1;
    size_t max_src_size_length = 0;
    size_t max_text_length = 0;
    for (auto& row : rows) {
        max_src_length = std::max(max_src_length, row.path_name.size());
        max_start_length = std::max(max_start_length, std::to_string(row.record_start).size());
        max_seq_size_length = std::max(max_seq_size_length, std::to_string(row.seq_size).size());
        max_src_size_length = std::max(max_src_size_length, std::to_string(row.path_length).size());
        max_text_length = std::max(max_text_length, row.aligned_seq.size());
    }
    // write and pad them
    for (auto& row : rows) {
        out << "s "
            << row.path_name << std::string(max_src_length - row.path_name.size(), ' ')
            << std::setw(max_start_length+1) << row.record_start
            << std::setw(max_seq_size_length+1) << row.seq_size
            << std::setw(max_is_rev_length+1) << (row.is_reversed ? " - " : " + ")
            << std::setw(max_src_size_length+1) << row.path_length
            << row.aligned_seq
            << "\n";
    }
    out << std::endl;
}
