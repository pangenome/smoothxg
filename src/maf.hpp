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
    std::pair<std::string, maf_partial_row_t> consensus_row;
};

void write_maf_row(std::ofstream &out, const maf_row_t& row){
    out << "s " + row.path_name + " " +
           std::to_string(row.record_start) + " " +
           std::to_string(row.seq_size) +
           (row.is_reversed ? " - " : " + ") +
           std::to_string(row.path_length) + " " +
           row.aligned_seq + "\n";
}
