#ifndef SMOOTHXG_RKMH_HPP
#define SMOOTHXG_RKMH_HPP

#include "mkmh.hpp"

// From Eric's https://github.com/edawson/rkmh

void hash_sequences(std::vector<std::string*> &seqs,
                    std::vector<std::vector<mkmh::hash_t>> &hashes,
                    std::vector<int> &hash_lengths,
                    std::vector<int> &kmer);

double compare(std::vector<mkmh::hash_t> alpha, std::vector<mkmh::hash_t> beta, int kmerSize);

#endif //SMOOTHXG_RKMH_HPP
