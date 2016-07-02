#ifndef EQUIV_D
#define EQUIV_D

#include <string>
#include <vector>
#include <tuple>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <omp.h>
#include "mkmh.hpp"


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
/*
 * int64_t [] make_hash_to_depth(min_depth, vector<int64_t> hashes);
 * */



inline string classify(vector<string>& read_kmers, map<string, vector<string> >& kmer_to_samples){

    
};

inline tuple<string, int, int> classify_and_count(vector<int64_t>& read_hashes, map<string, vector<int64_t> >& ref_to_hashes){
     //Parallel compare through the map would be really nice...
    int max_shared = 0;
    string sample = "";
    int shared_intersection = 0;
    int total_union = 0;
    map<string, vector<int64_t> >::iterator iter;
    for (iter = ref_to_hashes.begin(); iter != ref_to_hashes.end(); iter++){
         vector<int64_t> matches = hash_intersection(read_hashes, iter->second);
         //cerr << "MATCHES: " << matches.size() << endl;
         if (matches.size() > max_shared){
            max_shared = matches.size();
            sample = iter->first;
            shared_intersection = matches.size();
            total_union = read_hashes.size(); //hash_union(read_hashes, iter->second).size();

            //cerr << "Matches now: " << matches.size() << " " << sample << endl;
         }
    }
    return std::make_tuple(sample, shared_intersection, total_union);   
};


inline tuple<string, int, int> p_classify_and_count(vector<int64_t>& read_hashes, map<string, vector<int64_t> >& ref_to_hashes){
     //Parallel compare through the map would be really nice...
    int max_shared = 0;
    string sample = "";
    int shared_intersection = 0;
    int total_union = 0;
    vector<pair<string, vector<int64_t> > > ref_pairs(ref_to_hashes.begin(), ref_to_hashes.end());
    #pragma omp parallel for
    for (int i = 0; i < ref_pairs.size(); i++){
         vector<int64_t> matches = hash_intersection(read_hashes, ref_pairs[i].second);
         #pragma omp critical
         {
         if (matches.size() > max_shared){
            max_shared = matches.size();
            sample = ref_pairs[i].first;
            shared_intersection = matches.size();
            total_union = read_hashes.size(); //hash_union(read_hashes, iter->second).size();
         }
         }
    }
    return std::make_tuple(sample, shared_intersection, total_union);   
};


inline tuple<string, int, int> kmer_classify(vector<string>& readmers, map<string, vector<string> >& ref_mers){
    int max_shared = 0;
    string sample = "";
    int inter = 0;
    int uni = 0;
    map<string, vector<string> >::iterator iter;
    for (iter = ref_mers.begin(); iter != ref_mers.end(); iter++){
        vector<string> matches = kmer_intersection(readmers, iter->second);
        //cerr << "Matches: " << matches.size() << endl;
        if (matches.size() >= max_shared){
            max_shared = matches.size();
            inter = matches.size();
            sample = iter->first;
            uni = readmers.size();
        }
    }
     return std::make_tuple(sample, inter, uni);
};

inline string classify(vector<int64_t>& read_hashes, map<string, vector<int64_t> >& ref_to_hashes){

    //Parallel compare through the map would be really nice...
    int max_shared = 0;
    string ret = "";
    map<string, vector<int64_t> >::iterator iter;
    for (iter = ref_to_hashes.begin(); iter != ref_to_hashes.end(); iter++){
         vector<int64_t> matches = hash_intersection(read_hashes, iter->second);
         if (matches.size() > max_shared){
            ret = iter->first;
            max_shared = matches.size();
         }
    }
    return ret;
};

#endif
