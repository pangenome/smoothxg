#ifndef MKMH_D
#define MKMH_D

#include <vector>
#include <queue>
#include <set>
#include <unordered_map>
#include <string>
#include <cstring>
#include <sstream>
#include <locale>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <omp.h>
#include <assert.h>
#include <bitset>

#include "murmur3/murmur3.hpp"
#include "HASHTCounter.hpp"

#define DBGG

namespace mkmh {
    using namespace std;
    using namespace mkmh;

    typedef uint64_t hash_t;

    const static uint64_t MOD_HASH_ONE = 654;
    const static uint64_t MOD_HASH_TWO = 3459922;
    const static uint64_t MOD_HASH_THREE = 42;

    struct mkmh_kmer_list_t {
        char **kmers;
        int length;
        int k;

        mkmh_kmer_list_t() {

        };

        mkmh_kmer_list_t(int length, int k) {
            length = length;
            k = k;
            kmers = new char *[length];
        };

        ~mkmh_kmer_list_t() {
            for (int i = 0; i < length; ++i) {
                delete[] kmers[i];
            }
            delete[] kmers;
        };
    };

    struct mkmh_minimizer {
        uint64_t pos;
        uint32_t length;
        string seq;

        bool operator<(const mkmh_minimizer &rhs) const { return seq < rhs.seq; };
    };

    struct mkmh_hash_vec {
        hash_t *hashes;
        uint32_t size;
        uint64_t capacity;

        mkmh_hash_vec() {
            size = 0;
            capacity = 1000;
            hashes = new hash_t[capacity];
        };

        mkmh_hash_vec(uint32_t cap) {
            size = 0;
            capacity = cap;
            hashes = new hash_t[capacity];
        };

        mkmh_hash_vec(const mkmh_hash_vec &other) {
            size = other.size;
            capacity = other.capacity;
            hashes = new hash_t[capacity];
            for (size_t i = 0; i < other.size; ++i) {
                hashes[i] = other.hashes[i];
            }
        };

        ~mkmh_hash_vec() {
            if (this->capacity > 0) {
                delete[] this->hashes;
            }
        };

        inline void set_capacity(int cap) {
            if (this->capacity > 0) {
                delete[] hashes;
            }
            size = 0;
            capacity = cap;
            hashes = new hash_t[capacity];
        };

        void resize(double factor = 1.2) {
            size_t newcap = int(factor * capacity);
            hash_t *new_hashes = new hash_t[newcap];
            capacity = newcap;
            for (int i = 0; i < size; ++i) {
                new_hashes[i] = hashes[i];
            }
            delete[] hashes;
            hashes = new_hashes;
        };

        void emplace(hash_t h) {
            if (this->size == this->capacity) {
                resize();
            }
            this->hashes[this->size++] = h;
        };

        // Trim excess array space IFF the number of elements is
        // less than some percentage of capacity
        void trim(double factor = 0.6) {
            if (this->size < int(this->capacity * factor)) {
                int newcap = size;
                hash_t *new_hashes = new hash_t[newcap];
                for (int i = 0; i < size; ++i) {
                    new_hashes[i] = hashes[i];
                }
                delete[] hashes;
                hashes = new_hashes;
                capacity = newcap;
            }
        };

        void sort() {
            std::sort(this->hashes, this->hashes + this->size);
        };
    };

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

    static const int DNA_base_index[127] = {
            4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 0, 4, 1, 4, 4, 4,
            2, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 3, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 0, 4, 1, 4,
            4, 4, 2, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 3, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4
    };

    // Reverse complement lookup table
    static char rev_arr[26] = {
            84, 66, 71, 68, 69,
            70, 67, 72, 73, 74,
            75, 76, 77, 78, 79,
            80, 81, 82, 83, 65,
            85, 86, 87, 88, 89, 90
    };


    // Check a string (as a char*) for non-canonical DNA bases
    inline bool canonical(const char *x, int len) {
        bool trip = false;
        for (int i = 0; i < len; ++i) {
            trip |= valid_dna[x[i]];
        }
        return !trip;
    };

    inline bool canonical(string seq) {
        const char *x = seq.c_str();
        int len = seq.length();
        return canonical(x, len);
    };

    /* Reverse complement the string seq
     * (assumes seq is DNA, and returns non-ACTG letters as-is*/

    /* Reverse complement a C string
     * NB: does not check safety of string lengths.
     * NB: ret is modified to hold the reverse complement of seq.
     * */

    inline void reverse_complement(const char *seq, char *ret, int len) {

        //assert(seq != ret);
        if (ret == NULL) {
            ret = new char[len + 1];
        }

        for (int i = len - 1; i >= 0; i--) {
            ret[len - 1 - i] = (char) rev_arr[(int) seq[i] - 65];
        }
        ret[len] = '\0';
    };

    /* Reverse complement a string */
    inline string reverse_complement(string &seq) {
        const char *s = seq.c_str();
        int seqlen = seq.length();
        char *ret = new char[seqlen];
        reverse_complement(s, ret, seqlen);
        string s_revc(ret);
        delete[] ret;

        return s_revc;
    };

    /* Capitalize all characters in a string */
    /* Capitalize a C string */
    inline void to_upper(char *&seq, size_t length) {
        for (int i = 0; i < length; i++) {
            char c = seq[i];
            seq[i] = ((c - 91) > 0 ? c - 32 : c);
        }
    };

    inline void to_upper(const char *seq, int length, char *ret) {
        for (int i = 0; i < length; i++) {
            char c = seq[i];
            ret[i] = ((c - 91) > 0 ? c - 32 : c);
        }
    };

    /* Capitalize a string */
    inline string to_upper(string &seq) {
        for (int i = 0; i < seq.length(); i++) {
            char c = seq[i];
            seq[i] = ((c - 91) > 0 ? c - 32 : c);
        }
        return seq;
    };

    /* Reverse a string */
    inline string reverse(string seq) {
        string copy = string(seq);
        std::reverse(copy.begin(), copy.end());
        return copy;
    };

    /* Reverse a C string*/
    inline void reverse(char *seq, const int &len) {
        char tmp;
        for (int i = 0; i < len; ++i) {
            tmp = seq[len - i];
            seq[len - i] = seq[i];
            seq[i] = tmp;
        }
    };

    /** Custom sort function wrapping the STL implementation, mostly to allow descending sort. */
    inline void sort(hash_t *&hashes, int len, bool descending = false) {
        if (!descending) {
            std::sort(hashes, hashes + len);
        } else {
            std::sort(hashes, hashes + len, std::less<uint64_t>());

        }
    };

    inline void sort(vector <hash_t> &hashes, bool descending = false) {
        if (!descending) {
            std::sort(hashes.begin(), hashes.end());
        } else {
            std::sort(hashes.begin(), hashes.end(), std::less<uint64_t>());
        }
    };


    inline void kmerize(char *seq, const int &seq_len, const int &k, char **kmers, int &kmer_num) {
        char **ret = new char *[seq_len - k];
        kmer_num = seq_len - k;
        for (int i = 0; i < kmer_num; ++i) {
            ret[i] = new char[k + 1];
            memcpy(ret[i], seq + i, k);
        }
    }

    /* Returns the forward kmers of a sequence */
    inline vector <string> kmerize(string seq, int k) {
        vector <string> ret(seq.length() - k, "");

#pragma omp parallel for
        for (int i = 0; i < seq.length() - k; i++) {
            string s = seq.substr(i, k);
            //#pragma omp atomic read
            ret[i] = s;
            //ret.push_back(s);
            //ret.push_back(reverse(reverse_complement(s)));
        }
        return ret;
    };

    inline mkmh_kmer_list_t kmerize(char *seq, int seq_len, int k) {
        mkmh_kmer_list_t ret;
        ret.kmers = new char *[seq_len - k];
        ret.k = k;
        ret.length = seq_len - k;

        for (int i = 0; i < ret.length; ++i) {
            char *km = new char[k + 1];
            memcpy(km, seq + i, k);
            ret.kmers[i] = new char[k + 1];
            ret.kmers[i] = km;
            ret.kmers[i][k] = '\0';
        }
        return ret;
    };


    inline bool strcompare(const char *a, const int &alen, const char *b, const int &blen) {
        if (alen != blen) {
            return false;
        }

        for (int i = 0; i < alen; ++i) {
            if (a[i] != b[i]) {
                return false;
            }
        }

        return true;
    };

    inline void count_kmer_occurrence(const char *seq, const int &seq_len, const char *kmer, const int &k, int &count) {

        count = 0;
        int k_num = seq_len - k;
        for (int i = 0; i < k_num; ++i) {
            int addit = strcompare((const char *) seq + i, k, kmer, k) ? 1 : 0;
            count += addit;
        }
    };


    inline void
    count_substring_occurrence(const char *seq, const int &seq_len, const char *kmer, const int &k, int &count) {

        count = 0;
        int s_num = seq_len - k;
        int i = 0;
        while (i < s_num) {
            if (strcompare(seq + i, k, kmer, k)) {
                ++count;
                i = i + k;
            } else {
                i += 1;
            }
        }
    };


    /* Print the kmers of a string, tab separated, to cout 
     *   avoids allocating any new memory. */
    inline void print_kmers(char *seq, const int &len, int k) {
        int kmerized_length = len - k;
        stringstream st;
        for (int i = 0; i < kmerized_length - 1; ++i) {
            int j = 0;
            while (j < k) {
                st << seq[i + j];
                ++j;
            }
            st << "\t";
        }
        int j = 0;
        while (j < k) {
            st << seq[kmerized_length - 1 + j];
            ++j;
        }
        st << endl;
        cout << st.str();
        st.str("");
    };


    /* Returns the forward and reverse-reverse complement kmers for all kmer sizes in k */
    inline vector <string> multi_kmerize(string seq, vector<int> kSizes) {
        int i = 0;
        vector <string> ret;
        //ret.reserve(kSizes.size() * 1000);
        for (auto k : kSizes) {
            vector <string> kmers = kmerize(seq, k);
            ret.reserve(ret.size() + kmers.size());
            ret.insert(ret.end(), kmers.begin(), kmers.end());

            //for (i = 0; i + k < seq.length(); i++){
            //    ret.push_back(seq.substr(i, i+k));
            //    ret.push_back(reverse(reverse_complement(seq.substr(i, i+k))));
            //}


        }
        return ret;
    };

    /* Returns a deduplicated set of string kmers */
    inline vector <string> kmer_set(vector <string> kmers) {
        set <string> uniqs = set<string>(kmers.begin(), kmers.end());
        vector <string> ret = vector<string>(uniqs.begin(), uniqs.end());
        return ret;
    };

    /* Returns a heap (priority queue) of the kmers of a read converted to ints. */

    /* Returns a heap (priority queue) of the kmers of the read */
    inline priority_queue <string> kmer_heap(string seq, vector<int> kmer) {

        vector <string> base;

        //priority_queue<string> ret(base.begin(), base.end());
        for (auto k : kmer) {
            vector <string> outmers(seq.length() - k, "");
            for (int i = 0; i < seq.length() - k; i++) {
                string forward = seq.substr(i, k);
                string revrev = reverse(reverse_complement(forward));

                //ret.push( (revrev < forward ? revrev : forward) );
                outmers[i] = (revrev < forward ? revrev : forward);
            }
            base.reserve(outmers.size() + base.size());
            base.insert(base.end(), outmers.begin(), outmers.end());
        }

        priority_queue <string> ret(base.begin(), base.end());
        return ret;
    };

    const static int first_bits[5] = {0, 0, 1, 1, 0};
    const static int second_bits[5] = {1, 0, 0, 1, 0};

    /* Converts a string kmer to an integer representation */
    inline bool kmer_to_integer(const char *kmer, const int &length, hash_t *&h) {
        assert(length < 31);
        *h = 0;
        std::bitset<64> rb;
        int bit = -1;
        for (int i = 0; i < length; ++i) {
            int index = DNA_base_index[kmer[i]];
            if (index == 4) {
                *h = 0;
                return false;
            }
            ++bit;
            rb.set(63 - bit, first_bits[index]);
            ++bit;
            rb.set(63 - bit, second_bits[index]);
        }
        *h = rb.to_ullong();
        return true;
    }

    inline bool kmer_to_integer(const char *kmer, const int &length, hash_t &h) {
        assert(length < 31);
        h = 0;
        std::bitset<64> rb;
        int bit = -1;
        for (int i = 0; i < length; ++i) {
            int index = DNA_base_index[kmer[i]];
            if (index == 4) {
                h = 0;
                return false;
            }
            ++bit;
            rb.set(63 - bit, first_bits[index]);
            ++bit;
            rb.set(63 - bit, second_bits[index]);
        }
        h = rb.to_ullong();
        return true;
    }

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
    inline void calc_hashes(const char *seq, const int &len,
                            const int &k, hash_t *&hashes, int &numhashes) {
        char *reverse = new char[k + 1];
        uint32_t rhash[4];
        uint32_t fhash[4];
        //hash_t tmp_fwd;
        //hash_t tmp_rev;
        numhashes = len - k;
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


    /**
     * Thanks to https://stackoverflow.com/questions/39242932/how-to-encode-char-in-2-bits
     * NB: limited to kmers < k = 32
     */
    inline uint64_t kmer_to_integer(const char *kmer, const int &length) {
        uint64_t ret = 0;
        uint64_t retrev = 0;

        char *rev;
        reverse_complement(kmer, rev, length);

        for (int i = 0; i < length; ++i) {
            ret = (ret << 2) | ((kmer[i] >> 1) & 3);
            retrev = (retrev << 2) | ((rev[i] >> 1) & 3);
        }

        return ret < retrev ? ret : retrev;
    };

    /* Returns a deduplicated set of kmers or hashes as a vector<T> */
    template<typename T>
    inline vector <T> v_set(vector <T> kmers) {
        set <T> s = set<T>(kmers.begin(), kmers.end());
        vector <T> ret = vector<T>(s.begin(), s.end());
        return ret;
    }

    /* Returns only the forward shingles size k of a sequence */
    inline vector <string> shingle(string seq, int k) {
        int i = 0;
        vector <string> ret;
        for (i = 0; i < seq.length() - k; i++) {
            ret.push_back(seq.substr(i, k));
        }
        return ret;
    }

    /** Returns an mkmh_minimizer struct, equivalent to a tuple(kmer, position, kmer length), for every position in the genome **/
    inline vector <mkmh_minimizer> kmer_tuples(string seq, int k) {
        vector <string> kmers = kmerize(seq, k);
        vector <mkmh_minimizer> tups(kmers.size());
        for (int i = 0; i < kmers.size(); i++) {
            mkmh_minimizer mm;
            mm.seq = kmers[i];
            mm.pos = i;
            mm.length = k;
            tups[i] = mm;
        }

        return tups;
    }

    /** Finds the (w, k) minimizers of a string **/
    inline vector <mkmh_minimizer> minimizers(string seq, int k, int w) {
        vector <mkmh_minimizer> ret;
        vector <mkmh_minimizer> kmert = kmer_tuples(seq, k);
        int i = 0;
        for (i = 0; i + w < kmert.size(); ++i) {
            // get and sort kmers in window (i, i + w)
            vector <mkmh_minimizer> window_kmers(kmert.begin() + i, kmert.begin() + i + w);
            std::sort(window_kmers.begin(), window_kmers.end());
            // TODO filter minimizers if needed, e.g. to remove poly-As
            // create an mkmh_minimizer struct
            // tuck minimizer in ret
            ret.push_back(*(window_kmers.begin()));
        }
        return v_set(ret);
    };

    /** Finds the (w,k) minimizers and reports all of them (including duplicates) **/
    inline vector <mkmh_minimizer> unreduced_minimizers(string seq, int k, int w) {
        vector <mkmh_minimizer> ret;
        vector <mkmh_minimizer> kmert = kmer_tuples(seq, k);
        int i = 0;
        for (i = 0; i + w < kmert.size(); ++i) {
            // get and sort kmers in window (i, i + w)
            vector <mkmh_minimizer> window_kmers(kmert.begin() + i, kmert.begin() + i + w);
            std::sort(window_kmers.begin(), window_kmers.end());
            // TODO filter minimizers if needed, e.g. to remove poly-As
            // create an mkmh_minimizer struct
            // tuck miimizer in ret
            ret.push_back(*(window_kmers.begin()));
        }
        return ret;
    };

    inline void top_64_bits(uint32_t *&hash_holder, hash_t &ret) {
        ret = static_cast<uint64_t>(hash_holder[0]) << 32 | hash_holder[1];
    };


    // Calculate a single hash of a sequence (usually a kmer).
    // Takes in:
    //      seq: a char* (not affected)
    //      len: the length of seq (not affected)
    //      reverse: a char* of length(seq); modified to hold reverse_complement of (seq)
    //      forhash: a 128-bit int (e.g. uint32_t[4]) for holding forward hash.
    //      revhash: a 128-bit int (e.g. uint32_t[4]) for holding reverse hash.
    //      finhash: a 64-bit int (e.g. uint64_t) which holds the lesser of (forhash, revhash);
    inline void calc_hash(const char *seq,
                          const int &len,
                          char *&reverse,
                          uint32_t *&forhash,
                          uint32_t *&revhash,
                          hash_t *&fin_hash) {


        if (canonical(seq, len)) {
            reverse_complement(seq, reverse, len);
            MurmurHash3_x64_128(seq, len, 42, forhash);
            MurmurHash3_x64_128(reverse, len, 42, revhash);
            hash_t tmp_fwd = static_cast<uint64_t>(forhash[0]) << 32 | forhash[1];
            hash_t tmp_rev = static_cast<uint64_t>(revhash[0]) << 32 | revhash[1];
            *fin_hash = (tmp_fwd < tmp_rev ? tmp_fwd : tmp_rev);
        } else {
            *forhash = 0;
            *revhash = 0;
            *fin_hash = 0;
        }
    };


    /* Calculate the 64-bit hash for a string defined by seq and the length of seq */
    inline hash_t calc_hash(const char *seq, int seqlen) {
        char *reverse = new char[seqlen];
        uint32_t *fhash = new uint32_t[4];
        uint32_t *rhash = new uint32_t[4];
        hash_t *fin_hash = new hash_t[1];
        calc_hash(seq, seqlen, reverse, fhash, rhash, fin_hash);
        hash_t ret = *(fin_hash);
        delete[] reverse;
        delete[] fhash;
        delete[] rhash;
        delete[] fin_hash;
        return ret;
    }

    /* Calculate the hash for a string seq */
    inline hash_t calc_hash(string seq) {
        int k = seq.length();
        const char *x = seq.c_str();
        return calc_hash(x, k);
    };




    /** calc_hashes for multiple kmers sizes **/
    /** returns an array, with hashes for each kmer size concatenated
     * to those of the previous kmer size
     **/
    inline void calc_hashes(const char *seq, int seq_length,
                            vector<int> kmer_sizes,
                            hash_t *&hashes, int &numhashes) {
        numhashes = 0;

        // This holds the number of hashes preceeding the
        // kmer size currently being hashed.
        vector<int> offsets;
        for (auto k : kmer_sizes) {
            offsets.push_back(numhashes);
            numhashes += seq_length - k;
        }
        hashes = new hash_t[numhashes];

        for (int i = 0; i < kmer_sizes.size(); ++i) {
            int k = kmer_sizes[i];
            int local_numhash;
            hash_t *l_start;
            calc_hashes(seq, seq_length, k, l_start, local_numhash);
            memcpy(hashes + offsets[i], l_start, local_numhash * sizeof(hash_t));
            delete[] l_start;
        }
    }

    inline void calc_hashes(const char *seq, const int &len,
                            const int &k, hash_t *&hashes, int &numhashes, mkmh::HASHTCounter *&htc) {
        char *reverse = new char[k + 1];
        uint32_t rhash[4];
        uint32_t fhash[4];
        //hash_t tmp_fwd;
        //hash_t tmp_rev;
        numhashes = len - k;
        hashes = new hash_t[numhashes];
        for (int i = 0; i < numhashes; ++i) {
            if (canonical(seq + i, k)) {
                reverse_complement((seq + i), reverse, k);
                MurmurHash3_x64_128((seq + i), k, 42, fhash);
                MurmurHash3_x64_128(reverse, k, 42, rhash);
                //hash_t tmp_fwd = *((hash_t*) fhash);
                //hash_t tmp_rev = *((hash_t*) rhash);
                hash_t tmp_fwd = static_cast<uint64_t>(fhash[0]) << 32 | fhash[1];
                hash_t tmp_rev = static_cast<uint64_t>(rhash[0]) << 32 | rhash[1];

                hashes[i] = (tmp_fwd < tmp_rev ? tmp_fwd : tmp_rev);
                htc->increment(hashes[i]);
            } else {
                hashes[i] = 0;
            }

        }
        delete[] reverse;
    }

    inline void calc_hashes(const char *seq, int seq_length,
                            vector<int> kmer_sizes,
                            hash_t *&hashes, int &numhashes,
                            HASHTCounter *&htc) {

        numhashes = 0;

        // This holds the number of hashes preceeding the
        // kmer size currently being hashed.
        vector<int> offsets;
        for (auto k : kmer_sizes) {
            offsets.push_back(numhashes);
            numhashes += seq_length - k;
        }
        hashes = new hash_t[numhashes];

        for (int i = 0; i < kmer_sizes.size(); ++i) {
            int k = kmer_sizes[i];
            int local_numhash;
            //hash_t* l_start = hashes + offsets[i];
            hash_t *l_start;
            // HTC gets incremented within this function, so no need to do a bulk increment.
            calc_hashes(seq, seq_length, k, l_start, local_numhash, htc);
            memcpy(hashes + offsets[i], l_start, local_numhash * sizeof(hash_t));
            delete[] l_start;
        }
    };


    /* Calculate all the hashes of the kmers length k of seq */
    inline vector <hash_t> calc_hashes(const char *seq, int seq_length, int k) {
        int numhashes = 0;
        hash_t *hashes;
        calc_hashes(seq, seq_length, k, hashes, numhashes);
        vector <hash_t> ret(numhashes);
        for (int i = 0; i < numhashes; i++) {
            ret[i] = *(hashes + i);
        }

        delete[] hashes;

        return ret;
    };

    /* Calculate all the hashes of the kmers length k of seq */
    inline vector <hash_t> calc_hashes(string seq, int k) {
        const char *x = seq.c_str();
        int l = seq.length();
        return calc_hashes(x, l, k);
    }

    inline vector <hash_t> calc_hashes(const char *seq, const int &len, const vector<int> &k_sizes) {
        vector <hash_t> ret;
        for (auto k : k_sizes) {
            vector <hash_t> t = calc_hashes(seq, len, k);
            ret.insert(ret.end(), t.begin(), t.end());
        }
        return ret;
    };

    /** Calculate the hashes of seq
     *  and fill in a HASHTCounter htc so that
     *  hashes can be kept or removed based on the number of times
     *  they occur in seq.
     **/
    // void calc_hashes(const char* seq, const int& len,
    //         const int& k, hash_t*& hashes, int& numhashes, HASHTCounter*& htc);

    // void calc_hashes(const char* seq, const int& len,
    //         const int& k, hash_t*& hashes, int& numhashes, unordered_map<hash_t, int> counts);

    /** Calculate the hashes for kmers of multiple lengths in <kmer>
    */
    inline vector <hash_t> calc_hashes(string seq, const vector<int> &k_sizes) {
        const char *x = seq.c_str();
        int l = seq.length();
        return calc_hashes(x, l, k_sizes);
    };

    inline void
    hash_intersection_size(const hash_t *alpha, const int &alpha_size, const hash_t *beta, const int &beta_size,
                           int &ret) {
        int a_ind = 0;
        int b_ind = 0;
        ret = 0;
        while (a_ind < alpha_size && b_ind < beta_size) {
            if (alpha[a_ind] == beta[b_ind] && alpha[a_ind] != 0) {
                ++ret;
                ++a_ind;
                ++b_ind;
            } else if (alpha[a_ind] > beta[b_ind]) {
                ++b_ind;
            } else {
                ++a_ind;
            }
        }

    };

    inline void
    hash_set_intersection_size(const hash_t *alpha, const int &alpha_size, const hash_t *beta, const int &beta_size,
                               int &ret) {
        int a_ind = 0;
        int b_ind = 0;
        ret = 0;
        hash_t prev = 0;
        while (a_ind < alpha_size && b_ind < beta_size) {
            if (alpha[a_ind] == beta[b_ind] && alpha[a_ind] != prev && alpha[a_ind] != 0) {
                ++ret;
                prev = alpha[a_ind];
                ++a_ind;
                ++b_ind;

            } else if (alpha[a_ind] > beta[b_ind]) {
                ++b_ind;
            } else {
                ++a_ind;
            }
        }

    };

    /** Calculate a MinHash sketch for kmers length (2 * k) with skip bases in between the two k-length halves **/
    inline vector <hash_t> allhash_64_linkmer(string seq, int k, int skip) {
        vector <hash_t> ret(seq.size() - (k * 2));
        int last_kmer_ind = seq.size() - (skip + (2 * k));

#pragma omp parallel for
        for (int i = 0; i <= last_kmer_ind; ++i) {
            char *linkmer = new char[k * 2];
            for (int j = 0; j < k; ++j) {
                linkmer[j] = seq[i + j];
                linkmer[k + j] = seq[i + skip + k + j];
            }
            hash_t c_hash = calc_hash(linkmer, 2 * k);
            ret[i] = c_hash;
            delete[] linkmer;
        }

        return ret;
    }
    /* Returns the forward shingles of all k sizes of a sequence */
    /* Shingles are just forward-only kmers */
    inline vector <string> multi_shingle(string seq, vector<int> kSizes) {
        int i = 0;
        vector <string> ret;
        for (auto k : kSizes) {
            for (i = 0; i + k < seq.length(); i++) {
                ret.push_back(seq.substr(i, k));
            }
        }
        return ret;
    };

    /** Mask (by converting to zero) hashes that don't satisfy min_occ <= frequency(h) <= max_occ **/
    inline void mask_by_frequency(hash_t *&hashes, const int &num_hashes,
                                  HASHTCounter *htc,
                                  int min_occ = 0,
                                  uint32_t max_occ = UINT32_MAX) {

        for (int i = 0; i < num_hashes; ++i) {
            int freq = 0;
            htc->get(hashes[i], freq);
            hashes[i] = (min_occ <= freq && freq <= max_occ) ? hashes[i] : 0;
        }
    };

    inline vector <hash_t> minhashes(hash_t *hashes, int num_hashes, int sketch_size, bool useBottom) {
        vector <hash_t> x = vector<hash_t>(hashes, hashes + num_hashes);
        std::sort(x.begin(), x.end());

        int valid_ind = 0;
        while (x[valid_ind] == 0) {
            valid_ind++;
        }

        /*for (auto xx : x){
            if (xx != 0){
             cerr << xx << endl;
            }
        }*/

        int hashmax = valid_ind + sketch_size < num_hashes ? valid_ind + sketch_size : num_hashes - 1;
        return std::vector<hash_t>(x.begin() + valid_ind, x.begin() + hashmax);
    }

    /** MinHash - given an array of hashes, modify the mins array to hold 
     * the lowest/highest N (excluding zeros) **/
    inline void minhashes(hash_t *&hashes, int num_hashes,
                          int sketch_size,
                          hash_t *&ret,
                          int &retsize,
                          bool use_bottom = true) {

        int maxlen = min(num_hashes, sketch_size);
        ret = new hash_t[maxlen];
        retsize = 0;
        mkmh::sort(hashes, num_hashes, !use_bottom);

        int start = 0;
        while (retsize < sketch_size && start < num_hashes) {
            if (hashes[start] != 0) {
                ret[retsize] = hashes[start];
                ++retsize;
            }
            ++start;
        }
    };


    inline void minhashes_frequency_filter(hash_t *hashes, int num_hashes,
                                           int sketch_size,
                                           hash_t *&ret,
                                           int &retsize,
                                           HASHTCounter *htc,
                                           int min_occ = 0,
                                           uint32_t max_occ = UINT32_MAX,
                                           bool use_bottom = true) {

        ret = new hash_t[sketch_size];
        mkmh::sort(hashes, num_hashes, !use_bottom);

        int maxlen = min(num_hashes, sketch_size);
        int start = 0;
        while (retsize < sketch_size && start < num_hashes) {
            int freq = htc->get(hashes[start]);
            if (hashes[start] != 0 && freq <= max_occ && freq >= min_occ) {
                ret[retsize] = hashes[start];
                ++retsize;
            }
            ++start;
        }
    };

    inline void minhashes_min_occurrence_filter(hash_t *hashes, int num_hashes,
                                                int sketch_size,
                                                hash_t *&ret,
                                                int &retsize,
                                                HASHTCounter *htc,
                                                int min_occ = 0,
                                                bool use_bottom = true) {

        ret = new hash_t[sketch_size];
        mkmh::sort(hashes, num_hashes, !use_bottom);

        int maxlen = min(num_hashes, sketch_size);
        int start = 0;
        while (retsize < sketch_size && start < num_hashes) {
            if (hashes[start] != 0 && htc->get(hashes[start]) >= min_occ) {
                ret[retsize] = hashes[start];
                ++retsize;
            }
            ++start;
        }
    };

    /** Base MinHash function - return the lowest n = min(num_hashes, sketch_size) hashes. **/
    vector <hash_t> minhashes(hash_t *hashes, int num_hashes, int sketch_size, bool useBottom = true);

    /* Returns the lowest hashSize hashes of the kmers (length k...k` in k) of seq */
    inline vector <hash_t> minhash_64(string &seq, vector<int> &k, int hashSize, bool useBottom) {
        vector <hash_t> ret;
        //ret.reserve(k.size() * seq.size());

        for (auto km_sz : k) {
            vector <hash_t> tmp = calc_hashes(seq, km_sz);
            ret.insert(ret.end(), tmp.begin(), tmp.end());
        }
        std::sort(ret.begin(), ret.end());
        int nonzero_ind = 0;
        while (nonzero_ind < ret.size() && ret[nonzero_ind] == 0){
            nonzero_ind++;
        }

        hashSize += nonzero_ind;
        int hashmax = hashSize < ret.size() ? hashSize : ret.size() - 1;


        return useBottom ?
               vector<hash_t>(ret.begin() + nonzero_ind, ret.begin() + hashmax) :
               vector<hash_t>(ret.rbegin() + nonzero_ind, ret.rbegin() + hashmax);

    };

    /* Returns the bottom/top hashSize hashes of kmers size k in seq */
    inline vector <hash_t> minhash_64(string seq, int k, int hashSize, bool useBottom) {
        vector <string> kmers = kmerize(seq, k);
        vector <hash_t> ret(kmers.size(), 0);

        for (int i = 0; i < kmers.size(); i++) {

            uint32_t seed = 42;
            uint32_t khash[4];
            uint32_t rev_rev_khash[4];
            int str_length = kmers[i].size();
            const char *
                    forward = kmers[i].c_str();
            string rrf = reverse(reverse_complement(kmers[i]));
            const char *
                    rev_rev_forward = rrf.c_str();

            MurmurHash3_x64_128(forward, str_length, seed, khash);
            MurmurHash3_x64_128(rev_rev_forward, str_length, seed, rev_rev_khash);

            hash_t tmp_rev = *((hash_t *) rev_rev_khash);
            hash_t tmp_for = *((hash_t *) khash);

            ret[i] = tmp_for < tmp_rev ? tmp_for : tmp_rev;
        }

        std::sort(ret.begin(), ret.end());


        if (useBottom) {
            return vector<hash_t>(ret.begin(), ret.begin() + hashSize);
        } else {
            return vector<hash_t>(ret.rbegin(), ret.rbegin() + hashSize);
        }
    }

    /* helper function: returns the top hashSize hashes of the kmers size k in seq */
    inline vector <hash_t> top_minhash_64(string seq, int k, int hashSize) {
        return minhash_64(seq, k, hashSize, false);
    };


    /* helper function: returns the bottom hashSize hashes of the kmers size k in seq */
    inline vector <hash_t> bottom_minhash_64(string seq, int k, int hashSize) {
        return minhash_64(seq, k, hashSize, true);
    };

    /* Returns the bottom/top hashSize hashes of kmers size k which
     * occur more than minDepth times, based on the depth in hash_to_depth */
    inline vector <hash_t> minhash_64_depth_filter(string &seq, vector<int> &k,
                                                   int hashSize, bool useBottom, int minDepth,
                                                   unordered_map<hash_t, int> &hash_to_depth) {

        vector <string> kmers = multi_kmerize(seq, k);

        vector <hash_t> ret;
        ret.reserve(kmers.size());

        //#pragma omp parallel for
        for (int i = 0; i < kmers.size(); i++) {
            uint32_t khash[4];
            uint32_t rev_rev_khash[4];
            int str_length = kmers[i].size();
            const char *forward = kmers[i].c_str();

            string rrf = reverse(reverse_complement(kmers[i]));
            const char *rev_rev_forward = rrf.c_str();

            MurmurHash3_x64_128(forward, str_length, 42, khash);
            MurmurHash3_x64_128(rev_rev_forward, str_length, 42, rev_rev_khash);

            hash_t tmp_for = hash_t(khash[2]) << 32 | hash_t(khash[1]);
            hash_t tmp_rev = hash_t(rev_rev_khash[2]) << 32 | hash_t(rev_rev_khash[1]);

            hash_t r_hash = tmp_for < tmp_rev ? tmp_for : tmp_rev; //ret.push_back(r_hash);
            if (hash_to_depth[r_hash] > minDepth) {
#pragma omp critical
                ret.push_back(r_hash);
            }
        }

        if (ret.size() == 0) {
            return ret;
        }
        std::sort(ret.begin(), ret.end());

        int hashmax = hashSize < ret.size() ? hashSize : ret.size() - 1;

        return useBottom ?
               vector<hash_t>(ret.begin(), ret.begin() + hashmax) :
               vector<hash_t>(ret.rbegin(), ret.rbegin() + hashmax);
    };

    /* Takes in a list of pre-computed hashes and returns the MinHash (size hashSize)
     * of the hashes that pass the depth filter */

    inline vector <hash_t> minhash_64_depth_filter(vector <hash_t> &hashes, int hashSize, bool useBottom,
                                                   int min_depth, unordered_map<hash_t, int> &hash_to_depth) {
        vector <hash_t> ret;
        ret.reserve(hashes.size() / 2);
        for (int i = 0; i < hashes.size(); i++) {
            if (hash_to_depth[hashes[i]] > min_depth) {
#pragma omp critical
                ret.push_back(hashes[i]);
            }
        }

        /**
         * Special case if no hashes pass the depth filter.
         */
        if (ret.size() == 0) {
            return ret;
        }

        std::sort(ret.begin(), ret.end());

        int hashmax = hashSize < ret.size() ? hashSize : ret.size() - 1;
        return useBottom ?
               vector<hash_t>(ret.begin(), ret.begin() + hashmax) :
               vector<hash_t>(ret.rbegin(), ret.rbegin() + hashmax);
    };

    vector <hash_t> minhash_64_fast(string seq, vector<int> kmer, int sketchSize, bool isBottom = true);


    // inline void modimizer(char*& seq,
    //     const std::size_t& seqlen,
    //     const std::uint16_t& k,
    //     mkmh::mkmh_hash_vec*& ret,
    //     bool hashKmers = true){
    //     ret = new mkmh::mkmh_hash_vec(seqlen);
    // }

    /** TODO: this definitely has an off by W bug or something **/
    inline void minimizers(char *&seq,
                           const int &seqlen,
                           const int &k,
                           const int &w,
                           mkmh_hash_vec *&hvec,
                           bool hashKmers = true
    ) {
        hvec = new mkmh_hash_vec(((seqlen - k) / w) + ((seqlen - k) % w == 0 ? 0 : 1));
        vector<int> ksz;
        ksz.push_back(k);

        mkmh::hash_t *vals;
        int num_vals;

        calc_hashes(seq, seqlen, ksz, vals, num_vals);
        for (int i = 0; i < num_vals - w; i += w) {
            std::sort(vals + i, vals + i + w);
            hvec->emplace(vals[i]);
        }
    };

    /* Returns the intersection of alpha and beta, including duplicates */
    inline std::tuple<hash_t *, int> hash_intersection(hash_t alpha[], int alpha_start, int alpha_len,
                                                       hash_t beta[], int beta_start, int beta_len,
                                                       int sketch_size) {
        hash_t *ret = new hash_t[sketch_size];
        int ret_len = 0;

        int i = alpha_start;
        int j = beta_start;
        while (i < alpha_len && j < beta_len) {
            if (alpha[i] == beta[j]) {
                ret[ret_len] = alpha[i];
                ++i;
                ++j;
                ++ret_len;
            } else if (alpha[i] > beta[j]) {
                ++j;
            } else {
                ++i;
            }

        }

        //vector<hash_t> ret;
        //set_intersection(alpha + alpha_start, alpha + alpha_len,
        //                beta + beta_start, beta + beta_len,
        //                ret.begin());


        //return std::make_tuple(&(*ret.begin()), ret.size());
        return std::make_tuple(ret, ret_len);
    };

    /* Returns the intersection of alpha and beta, including duplicates the number of
       times they appear in both vectors */
    inline vector <hash_t> hash_intersection(vector <hash_t> alpha, vector <hash_t> beta) {
        vector <hash_t> ret;
        ret.reserve(alpha.size());
        int i = 0;
        int j = 0;
        while (i < alpha.size() && alpha[i] == 0){
            i++;
        }
        while(j < beta.size() && beta[j] == 0){
            j++;
        }
        while (i < alpha.size() && j < beta.size()) {
            if (alpha[i] == beta[j]) {
                ret.push_back(alpha[i]);
                i++;
                j++;
            } else if (alpha[i] > beta[j]) {
                j++;
            } else {
                i++;
            }
        }

        return ret;
    };

    /* Returns the union of the two sets after deduplicating all duplicates */
    inline vector <hash_t> hash_union(vector <hash_t> alpha, vector <hash_t> beta) {
        vector <hash_t> ret;


        ret.reserve(alpha.size() + beta.size());
        ret = vector<hash_t>(alpha.begin(), alpha.end());
        ret.insert(ret.end(), beta.begin(), beta.end());
        return ret;
    };

    /* Returns the union of the hashes in alpha and beta, including duplicates */
    inline vector <hash_t> hash_set_union(vector <hash_t> alpha, vector <hash_t> beta) {
        return v_set(hash_union(alpha, beta));
    };

    /* Returns the intersection of both sets. Duplicates are included only once */
    inline vector <hash_t> hash_set_intersection(vector <hash_t> alpha, vector <hash_t> beta) {
        return hash_intersection(v_set(alpha), v_set(beta));

    }

    inline void percent_identity(const hash_t *alpha, const int len,
                                 const hash_t *ref, const int reflen, double &ret) {
        int shared = 0;
        hash_intersection_size(alpha, len, ref, reflen, shared);
        ret = (double) shared / (double) len;
    };

    inline void sort_by_similarity(const hash_t *alpha, const int len,
                                   const vector <string> &refnames, int numrefs,
                                   const vector<hash_t *> &refhashes, const vector<int> &reflens,
                                   vector <string> &ret_names, vector<double> &ret_sims) {

        ret_names.resize(numrefs);
        ret_sims.resize(numrefs);

        vector <pair<int, double>> helper_vec(numrefs);
        for (int i = 0; i < numrefs; ++i) {
            double r = 0.0;
            percent_identity(alpha, len, refhashes[i], reflens[i], r);
            helper_vec[i] = std::make_pair(i, r);
        }
        sort(helper_vec.begin(), helper_vec.end(), [](const pair<int, double> &lhs, const pair<int, double> &rhs) {
            return lhs.second > rhs.second;
        });

        for (int i = 0; i < numrefs; ++i) {
            ret_names[i] = refnames[helper_vec[i].first];
            ret_sims[i] = helper_vec[i].second;
        }

    };

    inline void sort_by_similarity(const hash_t *alpha, const int len,
                                   const vector <string> &refnames, int numrefs,
                                   const vector<hash_t *> &refhashes, const vector<int> &reflens,
                                   vector <string> &ret_names, vector<double> &ret_sims,
                                   vector<int> &intersection_sizes) {

        ret_names.resize(numrefs);
        ret_sims.resize(numrefs);
        intersection_sizes.resize(numrefs);

        vector <pair<int, int>> helper_vec(numrefs);

        for (int i = 0; i < numrefs; ++i) {
            int shared = 0;
            hash_intersection_size(alpha, len, refhashes[i], reflens[i], shared);
            helper_vec[i] = std::make_pair(i, shared);
        }
        sort(helper_vec.begin(), helper_vec.end(), [](const pair<int, int> &lhs, const pair<int, int> &rhs) {
            return lhs.second > rhs.second;
        });

        for (int i = 0; i < numrefs; ++i) {
            ret_names[i] = refnames[helper_vec[i].first];
            intersection_sizes[i] = helper_vec[i].second;
            ret_sims[i] = (double) helper_vec[i].second / (double) len;

        }

    };

    /* Returns two vectors, one of sequence names and one of percent similarity, sorted by percent similarity to alpha.
     * NB: input vectors should be sorted. */
    inline tuple <vector<string>, vector<double>> sort_by_similarity(vector <hash_t> alpha,
                                                                     vector <vector<hash_t>> comps,
                                                                     vector <string> comp_names) {
        vector <string> ret_names;
        vector<double> ret_sims;
        vector <pair<int, double>> helper_vec;

        for (int i = 0; i < comps.size(); ++i) {
            int divisor = alpha.size();
            int shared = hash_intersection(alpha, comps[i]).size();
            double pct_id = (double) shared / (double) divisor;
            helper_vec.push_back(make_pair(i, pct_id));
        }
        sort(helper_vec.begin(), helper_vec.end(), [](const pair<int, double> &lhs, const pair<int, double> &rhs) {
            return lhs.second > rhs.second;
        });

        for (int i = 0; i < helper_vec.size(); i++) {
            ret_names.push_back(comp_names[helper_vec[i].first]);
            ret_sims.push_back(helper_vec[i].second);
        }

        return make_tuple(ret_names, ret_sims);
    }
    //inline void sort_by_similarity(hash_t* alpha, int hashnum, hash_t** comps, int* comp_nums, vector<char*> comp_names)
    /** Return the intersection of two kmer heaps **/
    inline priority_queue <string> kmer_heap_intersection(priority_queue <string> alpha, priority_queue <string> beta) {
        vector <string> base;
        base.reserve(alpha.size());
        while (alpha.size() != 0 && beta.size() != 0) {
            if (alpha.top() == beta.top()) {
                base.push_back(alpha.top());
                alpha.pop();
                beta.pop();
            } else if (alpha.top() > beta.top()) {
                alpha.pop();
            } else if (alpha.top() < beta.top()) {
                beta.pop();
            }
        }

        priority_queue <string> ret(base.begin(), base.end());
        return ret;
    };

    /** Return the intersection of two lists of kmers **/
    inline vector <string> kmer_intersection(vector <string> alpha, vector <string> beta) {
        vector <string> ret;
        ret.reserve(alpha.size());
        int i = 0;
        int j = 0;
        while (i < alpha.size() && j < beta.size()) {
            if (alpha[i] == beta[j]) {
                ret.push_back(alpha[i]);
                i++;
                j++;
            } else if (alpha[i] > beta[j]) {
                j++;
            } else {
                i++;
            }
        }

        return ret;
    };

    // Hold a buffer of bufsz kmers and pass that to the stringstream
    inline void print_kmers(char *seq, const int &len, int k, char *opt_char = NULL) {
        int kmerized_length = len - k;
        stringstream st;

        char *buf = new char[(k + 1)];
        buf[k] = '\0';

        if (opt_char != NULL) {
            st << opt_char << '\t';
        }
        for (int i = 0; i < kmerized_length - 1; ++i) {
            int j = 0;
            while (j < k) {
                buf[j] = seq[i + j];
                ++j;
            }
            st << buf << '\t';
        }
        int j = 0;
        while (j < k) {
            buf[j] = seq[kmerized_length - 1 + j];
            ++j;
        }
        st << buf;
        st << endl;
        cout << st.str();
        delete[] buf;
    };

    inline void print_hashes(hash_t *hashes, const int &numhashes, char *opt_char = NULL) {
        stringstream st;

        if (opt_char != NULL) {
            st << opt_char << '\t';
        }
        for (int i = 0; i < numhashes - 1; ++i) {
            st << hashes[i] << '\t';
        }
        st << hashes[numhashes - 1] << endl;
        cout << st.str();
    };

}

#endif
