#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include "murmur3.hpp"


namespace mkmh{
    using namespace std;

    /* Reverse the string seq */
    string reverse(string seq);
    
    /* Reverse complement the string seq (assumes seq is DNA, and returns non-ACTG letters as-is*/
    string reverse_complement(string seq);

    /* Returns the forward and reverse-reverse complement kmers of a sequence */
    vector<string> kmerize(string seq, int k);

    /* Returns the forward and reverse-reverse complement kmers for all kmer sizes in k */
    vector<string> multi_kmerize(string seq, vector<int> k);

    /* Returns a deduplicated set of string kmers */
    vector<string> kmer_set(vector<string> kmers);

    /* Returns a deduplicated set of kmers or hashes as a vector<T> */
    template<typename T>
    inline vector<T> v_set(vector<T> kmers){
        set<T> s = set<T>(kmers.begin(), kmers.end());
        vector<T> ret = vector<T>(s.begin(), s.end());
        return ret;
    }

    /* Returns the forward shingles size k of a sequence */
    vector<string> shingle(string seq, int k);

    /* Returns the forward shingles of all k sizes of a sequence */
    vector<string> multi_shingle(string seq, vector<int> k);

    /* Returns the lowest hashSize hashes of the kmers (length k...k` in k) of seq */
    vector<int64_t> minhash_64(string seq, vector<int> k, int hashSize, bool useBottom=true);

    /* Returns the bottom/top hashSize hashes of kmers size k in seq */ 
    vector<int64_t> minhash_64(string seq, int k, int hashSize, bool useBottom=true);

    /* helper function: returns the top hashSize hashes of the kmers size k in seq */
    vector<int64_t> top_minhash_64(string seq, int k, int hashSize);

    /* helper function: returns the bottom hashSize hashes of the kmers size k in seq */
    vector<int64_t> bottom_minhash_64(string seq, int k, int hashSize);

    /* Returns the union of the hashes in alpha and beta, including duplicates */
    vector<int64_t> hash_union(vector<int64_t> alpha, vector<int64_t> beta);

    /* Returns the intersection of alpha and beta, including duplicates the number of times they appear in both vectors */
    vector<int64_t> hash_intersection(vector<int64_t> alpha, vector<int64_t> beta);

    /* Returns the union of the two sets after deduplicating all duplicates */
    vector<int64_t> hash_set_union(vector<int64_t> alpha, vector<int64_t> beta);

    /* Returns the intersection of both sets. Duplicates are included only once */
    vector<int64_t> hash_set_intersection(vector<int64_t> alpha, vector<int64_t> beta);
}
