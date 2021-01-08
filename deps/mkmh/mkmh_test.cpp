#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "mkmh.hpp"
#include <fstream>
#include <string>
#include <algorithm>


using namespace std;
using namespace mkmh;

TEST_CASE("Reverse complement function works", "[reverse_complement]") {
    string t = "ACTGGCC";
    string rev = reverse_complement(t);

    SECTION("reverse_complement works on C++ strings.") {
        REQUIRE(rev == "GGCCAGT");
    }


    char k[6] = "AGGTC";
    char *ret = new char[6];
    char *retret = new char[6];
    reverse_complement(k, ret, 5);

    SECTION("reverse_complement returns expected string") {
        REQUIRE(strcmp(ret, "GACCT") == 0);
        REQUIRE(strlen(ret) == 5);
    }


    reverse_complement(ret, retret, 5);
    SECTION("reverse_complement does not affect its sequence pointer") {
        REQUIRE(ret == ret);
    }

    SECTION("reverse_complement, when applied twice, returns the original input string") {
        REQUIRE(*retret == *k);
    }


}

TEST_CASE("Canonical function works", "[canonical]") {
    string t = "ACTGGCNNNN";
    SECTION("canonical(string) catches Ns in a DNA string") {
        REQUIRE(canonical(t) == false);
    }

    char k[27] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    SECTION("canonical(char* , len) catches non-DNA letters") {
        REQUIRE(canonical(k, 26) == false);
    }

    char o[8] = "ACCCCTG";
    SECTION("canonical(char*, len) doesn't flag valid uppercase DNA characters") {
        REQUIRE(canonical(o, 7) == true);
    }


    char low[8] = "acccctg";
    SECTION("canonical(char*, len) doesn't flag valid lowercase DNA characters") {
        REQUIRE(canonical(low, 7) == true);
    }
}

TEST_CASE("Upper works for strings and chars") {
    char c[10] = "actgtgccc";
    char noncon[4] = "aBd";
    string d = "ABCDEFG";


}

TEST_CASE("v_set removes duplicates and returns a vector", "[v_set]") {

}

TEST_CASE("kmerize works as expected for strings", "[kmerize(string, ...)]") {

}

TEST_CASE("kmerize functions for char* work as expected", "[kmerize(char*, ...)]") {

}

TEST_CASE("minimizers behave as expected", "[minimizers]") {

}

TEST_CASE("Calc_hashes functions produce the right hashes", "[calc_hashes]") {

    char o[8] = "ACCCCTG";
    char t[8] = "ACCCCTG";

    SECTION("Hashes from calc_hashes for char* are consistent with those from calc_hash") {
        hash_t *h;
        int numhashes;
        calc_hashes((const char *) o, 7, 4, h, numhashes);

        for (int i = 0; i < numhashes; i++) {

            cerr << string(o + i, o + i + 4) << " : " << *(h + i) << " : " << calc_hash(o + i, 4) << endl;
        }

        bool trip = false;
        for (int i = 0; i < 7 - 4; i++) {
            trip = *(h + i) != calc_hash(o + i, 4);
        }
        REQUIRE(trip == false);
    }

    SECTION("calc_hashes functions all return the same hashes") {
        vector <hash_t> x = calc_hashes((const char *) o, 7, 4);
        hash_t *h;
        int num;

        calc_hashes((const char *) o, 7, 4, h, num);
        vector <hash_t> y(h, h + num);

        string zstr(o);
        vector <hash_t> z = calc_hashes(zstr, 4);

        vector<int> kmers;
        kmers.push_back(4);
        vector <hash_t> multis = calc_hashes(zstr, kmers);

        vector <hash_t> matched = calc_hashes((const char *) t, 7, 4);

        REQUIRE(std::mismatch(x.begin(), x.end(), y.begin(), y.end()).first == x.end());
        REQUIRE(std::mismatch(x.begin(), x.end(), z.begin(), z.end()).first == x.end());
        REQUIRE(std::mismatch(y.begin(), y.end(), z.begin(), z.end()).first == y.end());
        REQUIRE(std::mismatch(x.begin(), x.end(), multis.begin(), multis.end()).first == x.end());
        REQUIRE(std::mismatch(x.begin(), x.end(), matched.begin(), matched.end()).first == x.end());
    }

    SECTION("calc_hashes works with multiple kmer sizes") {
        vector<int> kmers;
        kmers.push_back(4);
    }

    SECTION("calc_hashes returns identical hashes for forward and reverse compliment") {
        char s[8] = "AAAAAAA";
        char t[8] = "TTTTTTT";
        vector <hash_t> s_hashes = calc_hashes(s, 7, 4);
        vector <hash_t> t_hashes = calc_hashes(t, 7, 4);
        REQUIRE(std::mismatch(s_hashes.begin(), s_hashes.end(), t_hashes.begin(), t_hashes.end()).first ==
                s_hashes.end());
    }

}

TEST_CASE("Sort and minhashes return the expect", "sort / minhash") {
    string x = "ACTGGCTTGCC";
    string y = "GGCAAGCCAGT";


}

TEST_CASE("Calc_hash family of functions work correctly", "[calc_hash()]") {
    string x = "ACTGGCTTGCC";
    string y = "GGCAAGCCAGT";

    SECTION("Hashes of forward and reverse-complement sequences are equal") {
        hash_t c_x = calc_hash(x);
        hash_t c_y = calc_hash(y);
        REQUIRE(c_x == c_y);

        REQUIRE(calc_hash("AAAAAA") == calc_hash("TTTTTT"));
    }

    SECTION("Hashes of calc_hash and calc_hashes are equivalent") {
        vector <hash_t> x_hashes = calc_hashes(x, 10);
        vector <hash_t> comp_hashes;
        for (int i = 0; i < x.length() - 10; i++) {
            comp_hashes.push_back(calc_hash(x.substr(i, i + 10)));
        }

        REQUIRE(std::mismatch(x_hashes.begin(), x_hashes.end(), comp_hashes.begin(), comp_hashes.end()).first ==
                x_hashes.end());
    }

    SECTION("Non-canonical bases cause a sequence to hash to zero") {
        string z = "ACGTNTTA";
        REQUIRE(calc_hash(z) == 0);
    }
}

TEST_CASE("hash_intersection family of functions work correctly", "[hash_intersection]") {

    hash_t *x = new hash_t[4];
    hash_t *y = new hash_t[6];
    int num;

    x[0] = 0;
    x[1] = 2;
    x[2] = 20938475420;
    x[3] = 987728;

    y[0] = 0;
    y[1] = 1;
    y[2] = 0;
    y[3] = 20938475420;
    y[4] = 10;
    y[5] = 987728;

    SECTION("fastest hash-intersection works") {
        hash_intersection_size(x, 4, y, 6, num);
        REQUIRE(num == 2);
    }

}

TEST_CASE("kmer_to_integer", "[kmer_to_integer]") {
    char s[8] = "ATAGAAA";
    char s_p[8] = "ATAGAAA";
    char non_s[13] = "ATAGAATTTTAA";
    char fail[12] = "ATAGANNNNAA";

    hash_t *z = new hash_t[1];
    hash_t *z_p = new hash_t[1];
    hash_t x;

    bool r = kmer_to_integer(s, 7, z[0]);
    kmer_to_integer(s_p, 7, z_p[0]);
    bool shouldfail = kmer_to_integer(fail, 11, x);
    hash_t nonz;

    kmer_to_integer(non_s, 12, nonz);
    REQUIRE(shouldfail == false);
    REQUIRE(x == 0);
    //cout << (hash_t) z << endl;
    REQUIRE(r == true);
    REQUIRE(z[0] == z_p[0]);
    REQUIRE(nonz != z[0]);

    delete[] z;
    delete[] z_p;
}
