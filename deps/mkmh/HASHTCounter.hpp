
#ifndef HTC_HPP
#define HTC_HPP

#include <cstdint>
#include <iostream>
#include <omp.h>
#include <cstdio>
#include <sstream>
#include <fstream>

namespace mkmh {
    using namespace mkmh;
    using namespace std;
    typedef uint64_t hash_t;

    class HASHTCounter {

    public:
        HASHTCounter();

        HASHTCounter(uint64_t sz);

        ~HASHTCounter();

        int &operator[](hash_t key);

        void increment(const hash_t &key);

        void bulk_increment(hash_t *h, int num);

        int &get(const hash_t &key);

        void get(const hash_t &key, int &ret);

        int size(void);

        void size(int sz);

        void resize(int sz);

        inline void set(int pos, int val) {
            *(counts + pos) = val;
        };

        int *begin(void);

        std::string to_string();

        void write_to_binary(std::string filename);

        void print();

    private:
        uint64_t my_size;
        int *counts;

    };
}

#endif
