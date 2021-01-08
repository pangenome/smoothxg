#include <vector>
#include <string>
#include "rkmh.hpp"
#include <cmath>
#include <algorithm>

void hash_sequences(std::vector<std::string *> &seqs,
                    std::vector<std::vector<mkmh::hash_t>> &hashes,
                    std::vector<int> &hash_lengths,
                    std::vector<int> &kmer) {
//#pragma omp parallel for
    for (int i = 0; i < seqs.size(); i++) {
        if (seqs[i] != nullptr) {
            hashes[i] = mkmh::calc_hashes(seqs[i]->c_str(), seqs[i]->length(), kmer);
            std::sort(hashes[i].begin(), hashes[i].end());
            hash_lengths[i] = hashes[i].size();
        }
    }
}

double compare(std::vector<mkmh::hash_t> alpha, std::vector<mkmh::hash_t> beta, int kmerSize) {
    int i = 0;
    int j = 0;

    uint64_t common = 0;
    uint64_t denom;

    while (i < alpha.size() && alpha[i] == 0) {
        i++;
    }
    while (j < beta.size() && beta[j] == 0) {
        j++;
    }
    denom = i + j;

    //todo early stopping
    while (i < alpha.size() && j < beta.size()) {
        if (alpha[i] == beta[j]) {
            i++;
            j++;
            common++;
        } else if (alpha[i] > beta[j]) {
            j++;
        } else {
            i++;
        }

        denom++;
    }

    // complete the union operation
    denom += alpha.size() - i;
    denom += beta.size() - j;

    //std::cerr << "common " << common << std::endl;
    //std::cerr << "denom " << denom << std::endl;

    double distance;

    //todo put a flag for denom: take the smallest between alpha.size, beta.size
    double jaccard = double(common) / denom;

    if (common == denom) // avoid -0
    {
        distance = 0;
    } else if (common == 0) // avoid inf
    {
        distance = 1.;
    } else {
        //distance = log(double(common + 1) / (denom + 1)) / log(1. / (denom + 1));
        distance = -log(2 * jaccard / (1. + jaccard)) / kmerSize;

        if (distance > 1) {
            distance = 1;
        }
    }

    return distance;
}
