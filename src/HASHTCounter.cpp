#include "HASHTCounter.hpp"

namespace HTC{ 
    using namespace HTC;
    using namespace std;

    HASHTCounter::HASHTCounter(){
        my_size = 1000000;
        counts = new int [my_size];
    }

    HASHTCounter::HASHTCounter(uint64_t sz){
        my_size = sz;
        counts = new int [my_size];
    }

    HASHTCounter::~HASHTCounter(){
        delete [] counts;
        my_size = 0;
    }

    string HASHTCounter::to_string(){
        stringstream sst;
        for (int i = 0; i < my_size; ++i){
            sst << counts[i] << endl;
        }
        return sst.str();
    }

    void HASHTCounter::print(){
        for (int i = 0; i < my_size; i++){
            cout << counts[i] << endl;
        }
    }

    void HASHTCounter::increment(htc_type key){
        //cout << (++counts [ key % my_size ]) << endl;
        #pragma omp atomic update
        ++(counts[ key % (uint64_t) my_size ]);
    }

    int& HASHTCounter::get(htc_type key){
        return (counts[ key % (uint64_t) my_size ]);
    }

    void HASHTCounter::get(htc_type key, int& ret){
        #pragma omp atomic write 
        ret = (counts[ key % (uint64_t) my_size ]);
    }

    int HASHTCounter::size(void){
        return my_size;
    }

    void HASHTCounter::size(int sz){
        
        delete [] counts;
        my_size = sz;
        counts = new int [my_size];

    }

    // TODO: not at all guaranteed safe.
    // Division / positioning in new array is uncheck, and wrong.
    void HASHTCounter::resize(int sz){
        
        int* n_counts = new int [ sz ];
        for (int i = 0; i < my_size; i++){
            *(n_counts + (i % sz)) = *(counts + i);    
        }
        my_size = sz;
        delete [] counts;
        counts = n_counts;
    }

    int& HASHTCounter::operator[](htc_type key){
        //value_t& operator[](std::size_t idx)       { return mVector[idx]; }
        //const value_t& operator[](std::size_t idx) const { return mVector[idx]; }
        return (counts [ key % my_size ]);
    }

    int* HASHTCounter::begin(void){
        return counts;
    }
}
