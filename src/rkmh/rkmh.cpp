#include <vector>
#include <string>
#include <iostream>
#include <unordered_set>
#include "rkmh.hpp"
#include "murmur3/murmur3.hpp"
#include <math.h>
#include <algorithm>

// Check a string (as a char*) for non-canonical DNA bases
inline bool canonical(const char *x, int len) {
    bool trip = false;
    for (int i = 0; i < len; ++i) {
        trip |= valid_dna[x[i]];
    }
    return !trip;
};
/* Reverse complement the string seq
 * (assumes seq is DNA, and returns non-ACTG letters as-is*/

/* Reverse complement a C string
 * NB: does not check safety of string lengths.
 * NB: ret is modified to hold the reverse complement of seq.
 * */

inline void reverse_complement(const char *seq, char *ret, int len) {

    //assert(seq != ret);
    if (ret == nullptr) {
        ret = new char[len + 1];
    }

    for (int i = len - 1; i >= 0; i--) {
        ret[len - 1 - i] = (char) rev_arr[(int) seq[i] - 65];
    }
    ret[len] = '\0';
};

/** Primary calc_hashes function **/
/** Takes the string to be hashed, its length,
 *  a single kmer size, a pointer to hold the hashes,
 *  and an integer to hold the number of hashes.
 *
 *  Possibly thread safe:
 *      seq, len and k are not modified
 *      new [] operator is known threadsafe
 *      User must handle hashes and numhashes properly in calling function.
 **/
inline void calc_hashes_(const char *seq, const uint64_t &len,
                         const uint64_t &k, hash_t *&hashes, int &numhashes) {
    char *reverse = new char[k + 1];
    uint32_t rhash[4];
    uint32_t fhash[4];
    //hash_t tmp_fwd;
    //hash_t tmp_rev;
    numhashes = (int) (len - k);
    hashes = new hash_t[numhashes];
    for (int i = 0; i < numhashes; ++i) {
        if (canonical(seq + i, k)) {
            reverse_complement(seq + i, reverse, k);
            MurmurHash3_x64_128(seq + i, k, 42, fhash);
            MurmurHash3_x64_128(reverse, k, 42, rhash);
            //hash_t tmp_fwd = *((hash_t*) fhash);
            //hash_t tmp_rev = *((hash_t*) rhash);
            hash_t tmp_fwd = static_cast<uint64_t>(fhash[0]) << 32 | fhash[1];
            hash_t tmp_rev = static_cast<uint64_t>(rhash[0]) << 32 | rhash[1];

            hashes[i] = (tmp_fwd < tmp_rev ? tmp_fwd : tmp_rev);
        } else {
            hashes[i] = 0;
        }
    }
    delete[] reverse;
};

/* Calculate all the hashes of the kmers length k of seq */
inline std::vector<hash_t> calc_hashes(const char *seq, uint64_t seq_length, uint64_t k) {
    int numhashes = 0;
    hash_t *hashes;
    calc_hashes_(seq, seq_length, k, hashes, numhashes);
    std::vector<hash_t> ret(numhashes);
    for (int i = 0; i < numhashes; i++) {
        ret[i] = *(hashes + i);
    }

    delete[] hashes;

    return ret;
};

inline std::vector<hash_t> calc_hashes(const char *seq, const uint64_t &len, const std::vector<uint64_t> &k_sizes) {
    std::vector<hash_t> ret;

    for (auto k : k_sizes) {
        //if (len > k) {
            std::vector<hash_t> t = calc_hashes(seq, len, k);
            ret.insert(ret.end(), t.begin(), t.end());
        //}
    }

    return ret;
};

void hash_sequences(std::vector<std::string *> &seqs,
                    std::vector<std::vector<hash_t>> &hashes,
                    std::vector<int> &hash_lengths,
                    std::vector<uint64_t> &kmer) {
//#pragma omp parallel for
    for (int i = 0; i < seqs.size(); i++) {
        if (seqs[i] != nullptr) {
            hashes[i] = calc_hashes(seqs[i]->c_str(), seqs[i]->length(), kmer);
            std::sort(hashes[i].begin(), hashes[i].end());
            hash_lengths[i] = hashes[i].size();
        }
    }
}

double compare(std::vector<hash_t> alpha, std::vector<hash_t> beta, uint64_t kmerSize) {
    int i = 0;
    int j = 0;

    uint64_t common = 0;
    uint64_t denom;

    while (alpha[i] == 0) {
        i++;
    }
    while (beta[j] == 0) {
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
