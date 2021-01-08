#include "HASHTCounter.hpp"

namespace mkmh {
    using namespace std;
    using namespace mkmh;

    HASHTCounter::HASHTCounter() {
        my_size = 1000000;
        counts = new int[my_size];
    }

    HASHTCounter::HASHTCounter(uint64_t sz) {
        my_size = sz;
        counts = new int[my_size];
    }

    HASHTCounter::~HASHTCounter() {
        delete[] counts;
        my_size = 0;
    }

    string HASHTCounter::to_string() {
        stringstream sst;
        for (int i = 0; i < my_size; ++i) {
            sst << counts[i] << endl;
        }
        return sst.str();
    }

    void HASHTCounter::write_to_binary(string filename) {
        ofstream ostr;
        ostr.open(filename, ios::out | ios::binary);
        ostr << this->my_size;
        for (int i = 0; i < my_size; ++i) {
            ostr << counts[i];
        }
        ostr.close();
    }

    void HASHTCounter::print() {
        for (int i = 0; i < my_size; i++) {
            cout << counts[i] << endl;
        }
    }

    void HASHTCounter::increment(const hash_t &key) {
        //cout << (++counts [ key % my_size ]) << endl;
#pragma omp atomic update
        ++(counts[key % static_cast<uint64_t>( my_size )]);
        /** #pragma omp critical
        {
             uint64_t k = key % (uint64_t) my_size;
             int v;
             v = *(counts + (int) k);
             v += 1;
             *(counts + (int) k)  = v;
         } */
    }

    void HASHTCounter::bulk_increment(hash_t *h, int num) {
        for (int i = 0; i < num; ++i) {
            this->increment(*(h + i));
        }
    }

    int &HASHTCounter::get(const hash_t &key) {
        return (counts[key % static_cast<uint64_t>(my_size)]);
    }

    void HASHTCounter::get(const hash_t &key, int &ret) {
        ret = (counts[key % static_cast<uint64_t>(my_size)]);
    }

    int HASHTCounter::size(void) {
        return my_size;
    }

    void HASHTCounter::size(int sz) {

        delete[] counts;
        my_size = sz;
        counts = new int[my_size];

    }

    // TODO: not at all guaranteed safe.
    // Division / positioning in new array is uncheck, and wrong.
    void HASHTCounter::resize(int sz) {

        int *n_counts = new int[sz];
        for (int i = 0; i < my_size; i++) {
            *(n_counts + (i % sz)) = *(counts + i);
        }
        my_size = sz;
        delete[] counts;
        counts = n_counts;
    }

    int &HASHTCounter::operator[](hash_t key) {
        //value_t& operator[](std::size_t idx)       { return mVector[idx]; }
        //const value_t& operator[](std::size_t idx) const { return mVector[idx]; }
        return (counts[key % my_size]);
    }

    int *HASHTCounter::begin(void) {
        return counts;
    }
}
