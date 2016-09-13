
#ifndef HTC_HPP
#define HTC_HPP
#include <cstdint>

namespace HTC{
    typedef uint64_t htc_type;
class HASHTCounter{

    public:
        HASHTCounter();
        HASHTCounter(int sz);
        ~HASHTCounter();
        int& operator[](htc_type key);

        void increment(htc_type key);
        int get(htc_type key);

        int size(void);
        void size(int sz);
        

    private:
        int my_size;
        int* counts;
        
};
}
#endif
