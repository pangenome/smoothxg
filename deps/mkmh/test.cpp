#include "mkmh.hpp"
#include <string>
#include <iostream>
#include <queue>

using namespace std;
using namespace mkmh;

bool testify(int t_num, string test, string obs, string exp) {
    if (obs == exp) {
        cout << "PASS: " << t_num << " " << test << endl;
        return true;
    } else {
        cout << "FAIL: " << t_num << " " << test << endl;
        return true;
    }
}

bool testify(int t_num, string test, bool x) {
    if (x) {
        testify(t_num, test, "", "");
        return true;
    } else {
        testify(t_num, test, "1", "2");
        return false;
    }
}

bool same(vector <hash_t> x, vector <hash_t> y) {

    bool ret = true;
    if (x.size() != y.size()) {
        return false;
    }

    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());
    for (int i = 0; i < x.size(); i++) {
        cerr << x[i] << "\t" << y[i] << endl;
        if (x[i] != y[i]) {
            ret = false;
            //return false;
        }
    }
    return ret;
    //return true;

}

int main() {
    int t_num = 0;

    string seq = "ATGCATGCATGCATGCATGC";
    vector <string> kmers = kmerize(seq, 5);
    bool x = true;
    for (auto e : kmers) {
        if (e.length() != 5) {
            x = false;
        }
    }

    /* 
     *
     * Test kmerize, shingle, and minhash64
     * 
     * */
    testify(t_num++, "The kmers produced by kmerize are the right size", x);

    vector <string> shingles = shingle(seq, 5);
    x = true;
    for (auto e : shingles) {
        if (e.length() != 5) {
            x = false;
        }
    }
    testify(t_num++, "The shingles produced by kmerize are the right size", x);

    vector <hash_t> ret = minhash_64(seq, 5, 5);
    testify(t_num++, "minhash_64 produces the right number of hashes", ret.size() == 5);
    vector <hash_t> o_ret = minhash_64(seq, 5, 5);
    testify(t_num++, "minhash_64 produces the same top and bottom values in the hash",
            (ret[0] == o_ret[0] & ret[4] == o_ret[4]));

    vector<int> ret_test = {5};
    ret = minhash_64(seq, ret_test, 5);
    testify(t_num++, "multi_kmer minhash_64 produces the right number of hashes", ret.size() == 5);
    for (int iii = 1; iii < 1000; iii++) {
        continue;
    }
    o_ret = minhash_64(seq, ret_test, 5);
    testify(t_num++, "multi_kmer minhash_64 produces the same top and bottom values in the hash",
            (ret[0] == o_ret[0] & ret[4] == o_ret[4]));


    vector <string> k_set = kmer_set(kmers);
    testify(t_num++, "Kmer set removes duplicate kmers", k_set.size() < kmers.size());

    /* Test top_minhash_64 and bottom_minhash_64 */
    x = true;
    vector <hash_t> tops = top_minhash_64(seq, 4, 4);
    vector <hash_t> bottoms = bottom_minhash_64(seq, 4, 4);
    for (auto e : tops) {
        for (auto f : bottoms) {
            if (e < f) {
                x = false;
            }
        }
    }
    for (int i = 0; i < bottoms.size(); i++) {
        //cerr << kmerize(seq, 4)[i] << endl;
        //cerr << bottoms[i] << " " << tops[i] << endl;
    }
    testify(t_num++, "top_minhash64 produces the bigger values than bottom_minhash_64", x);


    string seq2 = "ACTGaaatttt";
    vector<int> ks;
    ks.push_back(4);
    ks.push_back(4);
    kmers = kmerize(seq2, 4);
    vector <string> m_kmers = multi_kmerize(seq2, ks);

    /* Test multi_kmerize */
    testify(t_num++, "multi_kmerize produces twice as many kmers with two kmer sizes",
            kmers.size() * 2 == m_kmers.size());

    /* Test hash union / intersection */
    vector <hash_t> t1 = {1, 2, 3, 4, 5, 6};
    vector <hash_t> t2 = {4, 5, 6, 7, 8, 9};
    testify(t_num++, "Hash intersection of two sets is the expected size.", hash_intersection(t1, t2).size() == 3);

    t1 = {1, 2, 3, 4, 4, 4, 5, 6};
    t2 = {4, 4, 5, 6, 7, 8, 9};
    testify(t_num++, "Hash intersection counts duplicate values", hash_intersection(t1, t2).size() == 4);

    testify(t_num++, "Hash union yields the correct size with duplicate values",
            hash_union(t1, t1).size() == 2 * t1.size());
    //testify(t_num++ "Union / Intersection produces expected value", hash_intersection(t1, t1).size() / hash_intersection(
    testify(t_num++, "Hash union yields the correct size when duplicates are removed",
            v_set(hash_union(t1, t1)).size() == 6);

    testify(t_num++, "Hash set union produces the expected number of values", hash_set_union(t1, t2).size() == 9);
    testify(t_num++, "Hash set intersection produces the expected number of values",
            hash_set_intersection(t1, t2).size() == 3);

    string a = "AAATGCTTTTGCA";
    vector<int> three;
    three.push_back(3);
    priority_queue <string> a_heap = kmer_heap(a, three);
    //cerr << a_heap.top() << endl;
    testify(t_num++, "Kmer heap produces expected lowest kmer", (a_heap.top() == "CTT"));
    while (a_heap.size() > 1) {
        //cerr << a_heap.top() << endl;
        a_heap.pop();
    }
    testify(t_num++, "Kmer heap produces expected highest kmer", (a_heap.top() == "AAA"));

    //vector<hash_t> a_allhash_fast = allhash_unsorted_64_fast(a.c_str(), three);
    //testify(t_num++, "allhash_unsorted_64_fast produces a hash vector of the proper length", (a_allhash_fast.size() == 10));
    //for (auto ll : a_allhash_fast){
    //    cerr << ll << endl;
    //}

    vector<int> four;
    four.push_back(4);
    cerr << "Testing minhash 64 functions" << endl;
    vector <hash_t> bottoms_fast = minhash_64(seq, four, 4, true);
    testify(t_num++, "minhash_64 and minhash_64 produce the same hashes", same(bottoms, bottoms_fast));
    hash_t *fast_comp1 = new hash_t[4];

    fast_comp1[0] = 1234;
    fast_comp1[1] = 5534;
    fast_comp1[2] = 2312;
    fast_comp1[3] = 4532;
    hash_t *fast_comp2 = new hash_t[4];
    fast_comp2[0] = 1234;
    fast_comp2[1] = 2312;
    fast_comp2[2] = 1912;
    fast_comp2[3] = 3333;

    std::sort(fast_comp1, fast_comp1 + 4);
    std::sort(fast_comp2, fast_comp2 + 4);

    //std::tuple<hash_t*, int> hash_intersection(hash_t* alpha, int alpha_len, int alpha_start,
    //                                                                    hash_t* beta, int beta_len, int beta_start,
    //                               int sketch_size){

    std::tuple < hash_t * , int > inter = hash_intersection(fast_comp1, 0, 4, fast_comp2, 0, 4, 4);
    testify(t_num++, "new fast hash intersection produces the right number of values.", std::get<1>(inter) == 2);

    vector <mkmh_minimizer> mh = minimizers("ACTGGTTCCCCAAAATTTTTTGATAGTAGGATACCGACGC", 5, 10);
    testify(t_num++, "minimizers produces the right number of minimizers", mh.size() == 9);
    //cerr << "sz: " << mh.size() << endl;
    //for (auto z : mh){
    //cout << z.pos << " " << " " << z.length << " " << z.seq << endl;
    //}

    //cerr << seq << endl;
    vector <hash_t> link_hashes = allhash_64_linkmer(seq, 3, 0);
    vector<int> s;
    s.push_back(6);
    //vector<hash_t> base_kmers = allhash_unsorted_64(seq, s);
    //testify(t_num++, "linked_hashes produces the right number of hashes when skip = 0", link_hashes.size() == base_kmers.size());
    //bool tripped = false;
    //for (int i = 0; i < link_hashes.size(); i++){
    //    tripped = link_hashes[i] != base_kmers[i];        
    //}
    //testify(t_num++, "linked_hashes produces the correct hashes when skip = 0", !tripped);
    //char* eightseq = {'A', 'C', 'T', 'G', 'A', 'A', 'G', 'T', '\0'};
    //char eightseq [8] = {"ACTGAAG"};
    //print_kmers(eightseq, 7, 4);
    //cout << endl;

    //link_hashes = allhash_64_linkmer(seq, 4, 2);

    //char tt[] = "ATGGTTCCCGGTTTTT";
    //kmerize(tt, 16, 4);

    return 0;
}
