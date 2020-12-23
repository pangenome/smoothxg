#ifndef SMOOTHXG_RKMH_HPP
#define SMOOTHXG_RKMH_HPP

// From Eric's https://github.com/edawson/rkmh

#include <cstdint>

// Crazy hack char table to test for canonical bases
static const int valid_dna[127] = {
        1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 0, 1, 0, 1, 1, 1,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 0, 1, 0, 1,
        1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 0, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1
};

// Reverse complement lookup table
static char rev_arr[26] = {
        84, 66, 71, 68, 69,
        70, 67, 72, 73, 74,
        75, 76, 77, 78, 79,
        80, 81, 82, 83, 65,
        85, 86, 87, 88, 89, 90
};

typedef uint64_t hash_t;

void hash_sequences(std::vector<std::string> &keys,
                    std::vector<char *> &seqs,
                    std::vector<int> &lengths,
                    std::vector<std::vector<hash_t>> &hashes,
                    std::vector<int> &hash_lengths,
                    std::vector<int> &kmer);

std::vector<hash_t> hash_intersection(std::vector<hash_t> alpha, std::vector<hash_t> beta);

std::vector<hash_t> hash_union(std::vector<hash_t> alpha, std::vector<hash_t> beta);

#endif //SMOOTHXG_RKMH_HPP
