#ifndef EQUIV_D
#define EQUIV_D

#include <string>
#include "mkmh.hpp"
#include <vector>
#include <cstdint>
#include <map>

using namespace std;
using namespace mkmh;
inline map<string, vector<string> > make_kmer_to_samples(map<string, vector<string>>& sample_to_kmers){
    map<string, vector<string> > ret;
    map<string, vector<string> >::iterator sk_iter;
    for (sk_iter = sample_to_kmers.begin(); sk_iter != sample_to_kmers.end(); sk_iter++){
        for (auto km : sk_iter->second){
            ret[km].push_back(sk_iter->first);
        }
    }
    return ret;
};

inline map<string, vector<string> > make_sample_to_kmers(map<string, string>& name_to_sequence, int k){
    map<string, vector<string> > ret;

    map<string, string>::iterator ns_iter;
    for (ns_iter = name_to_sequence.begin(); ns_iter != name_to_sequence.end(); ns_iter++){
        ret[ns_iter->first] = kmerize(ns_iter->second, k);
    }


    return ret;

};

inline map<string, int> make_sample_to_count(vector<string>& read_kmers, map<string, vector<string> >& kmer_to_samples){
    map<string, int> sample_to_count;

    for (auto kmer : read_kmers){
        if (kmer_to_samples.count(kmer) != 0){
            vector<string> samples = kmer_to_samples[kmer];
            for (auto s : samples){
                sample_to_count[s] += 1;
            }
        }
    }

    return sample_to_count;
};

inline string classify(vector<string>& read_kmers, map<string, vector<string> >& kmer_to_samples){

    
};

struct Classification{
    string sample = "";
    int shared_intersection = 0;
    int total_union = 0;
};

inline struct Classification classify_and_count(vector<int64_t>& read_hashes, map<string, vector<int64_t> >& ref_to_hashes){
     //Parallel compare through the map would be really nice...
    int max_shared = -1;
    Classification ret;
    map<string, vector<int64_t> >::iterator iter;
    for (iter = ref_to_hashes.begin(); iter != ref_to_hashes.end(); iter++){
         vector<int64_t> matches = hash_set_intersection(read_hashes, iter->second);
         if (matches.size() > max_shared){
            ret.sample = iter->first;
            ret.shared_intersection = matches.size();
            ret.total_union = hash_set_union(read_hashes, iter->second).size();
         }
    }
    return ret;   
};

inline string classify(vector<int64_t>& read_hashes, map<string, vector<int64_t> >& ref_to_hashes){

    //Parallel compare through the map would be really nice...
    int max_shared = -1;
    string ret = "";
    map<string, vector<int64_t> >::iterator iter;
    for (iter = ref_to_hashes.begin(); iter != ref_to_hashes.end(); iter++){
         vector<int64_t> matches = hash_set_intersection(read_hashes, iter->second);
         if (matches.size() > max_shared){
            ret = iter->first;
            max_shared = matches.size();
         }
    }
    return ret;
};

#endif
