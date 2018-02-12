
#ifndef HTC_HPP
#define HTC_HPP
#include <cstdint>
#include <iostream>
#include <omp.h>
#include <cstdio>
#include <sstream>


namespace HTC{
    typedef uint64_t htc_type;
class HASHTCounter{

    public:
        HASHTCounter();
        HASHTCounter(uint64_t sz);
        ~HASHTCounter();
        int& operator[](htc_type key);

        void increment(htc_type key);
        int& get(htc_type key);

        void get(htc_type key, int& ret);

        int size(void);
        void size(int sz);
        void resize(int sz);
        inline void set(int pos, int val){
            *(counts + pos) = val;
        };

        int* begin(void);

        std::string to_string();
        void print();
        
    private:
        uint64_t my_size;
        int* counts;
        
};
}
#endif
