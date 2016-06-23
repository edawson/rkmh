# mkmh
Generate and compare MinHash signatures with multiple kmer sizes.

## Usage
To use mkmh functions in your code:  
1. Include the header file in your code  
    ```#include "mkmh.hpp"```      
2. Compile the library:  
    `` cd mkmh && make lib``  
3. Make sure the lib and header are on the LD include/lib paths (e.g. in your makefile):  
    `` gcc -o my_code my_code.cpp -L/path/to/mkmh -I/path/to/mkmh -lmkmh  
4. That's it!

## Available functions
The functions in the package produce vectors of hashes that are:  
1. Sorted  
2. Guaranteed 64 bits  
3. Bottom or Top hashes (i.e. the largest N hashes or smallest N hashes)  
4. From either kmers of size K or of all sizes [K<sub>1</sub> ... K<sub>m</sub>] provided.  


There are helper functions for ``top_minhash_64 `` and ``bottom_minhash_64``
as well as two `minhash_64` functions, one which takes a single kmer size and one which takes
a list of kmer sizes.


I tried to only use vectors in this library so there are functions for performing set operations on vectors.
This probably kills performance, but none of this has been tuned anyway so it's probably not
close to maximum efficiency.


## Getting help
Please reach out through (github)[https://github.com/edawson/mkmh] to post an issue,
shoot me an email, or file a bug report.
