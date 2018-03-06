#ifndef EQUIV_D
#define EQUIV_D

#include <string>
#include <set>
#include <vector>
#include <queue>
#include <tuple>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <omp.h>
#include "mkmh.hpp"
#include "HASHTCounter.hpp"


using namespace std;
using namespace mkmh;

struct c_comp{
    bool operator() (const hash_t& lhs, const hash_t& rhs) const {return lhs < rhs;};
};


struct min_heap{
    set<hash_t> con;
    int max_sz;
    void init(int s){
        max_sz = s;  
    };
    bool insert(hash_t val){
        if (val == 0){
            
        }
        else if (con.size() < 1){
           con.insert(val);
        }
        //else if (val < *(con.rbegin())){
        else{
           con.insert(val);
        }
    };
    void resize(){
        if (con.size() > max_sz){
            set<hash_t> x;
            int i = 0;
            auto itr = con.begin();
            while (i < max_sz){
                x.insert(*itr);
                i++;
                itr++;
            }
            con = x;
        }
    };
};


inline tuple<hash_t*, int*, int> merge(vector<vector<hash_t>> mins, int sketch_size, int max){
    HASHTCounter ht(100000000);
    min_heap k;
    k.init(sketch_size * 2);
    for (auto x : mins){
        for (auto y : x){
            ht.increment(y);
            k.insert(y);
        }   
    }
    k.resize();
    
    hash_t* ret = new hash_t [k.con.size()];
    int ret_sz = 0;
    int* counts = new int [k.con.size()];
    for (auto x = k.con.begin(); x != k.con.end(); x++){
       ret[ret_sz] = *x;
       counts[ret_sz] = ht.get(*x);
       ++ret_sz;
    }
    
    return std::make_tuple(ret, counts, ret_sz);
    

};

/**
 * takes a vector<vector<hasht>>, the mins of the reads set AND the sketch size
 * return a sketch of the lowest 
 */
inline tuple<hash_t*, int> merge(vector<vector<hash_t>> mins, int sketchsz){
    int ret_sz = 0;
    hash_t* ret = new hash_t [sketchsz];
    min_heap k;
    k.init(sketchsz);
    
    for (auto x : mins){
        //std::set<hash_t> y (x.begin(), x.end()); 
        for (auto y : x){
            k.insert(y);
        }
    }
    
    k.resize();
    for (auto x = k.con.begin(); x != k.con.end(); x++){
        ret[ret_sz++] = *x;
    }



    return make_tuple(ret, ret_sz);

};

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

/**
inline map<hash_t, int> pos_to_depth(string ref, vector<string> read_mers){

};
**/


/**
 vector<string> right_kmers(ref_kmers, read_kmers / sequence);
 * **/
inline map<hash_t, int> make_kmer_to_sample_count(map<string, hash_t*> name_to_hashes, map<string, int> name_to_num_hashes){
    map<hash_t, int> ret;
    map<hash_t, set<string> > helper;
    vector<pair<string, hash_t*> >  vp_name_hashes(name_to_hashes.begin(), name_to_hashes.end());
    #pragma omp for
    for (int i = 0; i < vp_name_hashes.size(); i++){
        for (int j = 0; j < name_to_num_hashes[vp_name_hashes[i].first]; j++){
            helper[vp_name_hashes[i].second[j]].insert(vp_name_hashes[i].first);
        }
    }

    for (auto x : helper){
        ret[x.first] = x.second.size();
    }


    return ret;

};


inline map<hash_t, int> make_kmer_to_sample_count(vector<pair<string, vector<hash_t> > > name_to_hashes){
    map<hash_t, int> ret;
    map<hash_t, set<string> > helper;
    for (int i = 0; i < name_to_hashes.size(); i++){
        for (int j = 0; j < name_to_hashes[i].second.size(); j++){
            helper[name_to_hashes[i].second[j]].insert(name_to_hashes[i].first);
        }
    }

    for (auto x : helper){
        ret[x.first] = x.second.size();
    }
    return ret;

};

inline map<string, vector<hash_t>> only_informative_kmers(map<string, hash_t*> name_to_hashes, map<string, int> name_to_num_hashes, int max_samples){
    map<hash_t, int> hash_to_sample_count = make_kmer_to_sample_count(name_to_hashes, name_to_num_hashes);

    map<string, vector<hash_t> > ret;
    for (auto x : name_to_hashes){
        for (int i = 0; i < name_to_num_hashes[x.first]; i++){
            if (hash_to_sample_count[x.second[i]] < max_samples){
                ret[x.first].push_back(x.second[i]);
            }
        }
    }

    return ret;

};

inline map<string, vector<hash_t> > only_informative_kmers(map<string, vector<hash_t> >& name_to_hashes, int max_samples){
    auto vv_count = vector<pair<string, vector<hash_t> > > (name_to_hashes.begin(), name_to_hashes.end());
    map<hash_t, int> hash_to_count = make_kmer_to_sample_count(vv_count);
    map<string, vector<hash_t> > ret;
    for (auto x : name_to_hashes){
        for (int i = 0; i < x.second.size(); i++){
            if (hash_to_count[x.second[i]] < max_samples){
                ret[x.first].push_back(x.second[i]);
            }
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


inline tuple<string, int, int> kmer_heap_classify(priority_queue<string> readmers, vector<pair<string, priority_queue<string> > > ref_mers){
    int max_shared = 0;
    string sample = "";
    int shared_intersection = 0;
    int total_union = 0;
    for (int i = 0; i < ref_mers.size(); i++){
        priority_queue<string> matches = kmer_heap_intersection(readmers, ref_mers[i].second);
        if (matches.size() > max_shared){
            max_shared = matches.size();
            sample = ref_mers[i].first;
            shared_intersection = matches.size();
            total_union = readmers.size();
        }
    }
    return std::make_tuple(sample, shared_intersection, total_union); 
};

inline tuple<string, int, int, bool> classify(vector<hash_t>& read_hashes, vector<int>& ref_counts, vector<pair<string, vector<hash_t> > >& ref_hashes, int min_diff){
    vector<int>::iterator max_iter = max_element(ref_counts.begin(), ref_counts.end());
    int max_val = int(*max_iter);
    int max_index = distance(ref_counts.begin(), max_iter);
    bool fail_min_diff = false;
    string sample = ref_hashes[max_index].first;
    sort(ref_counts.begin(), ref_counts.end(), std::greater<int>());
    if (ref_counts.size() >= 2 && ref_counts[0] - ref_counts[1] <= min_diff){
        fail_min_diff = true;
    }
    
    return make_tuple(sample, max_val, read_hashes.size(), fail_min_diff);

}

inline vector<int> all_count(vector<hash_t> read_hashes, vector<pair<string, vector<hash_t> > >& ref_to_hashes){
     //Parallel compare through the map would be really nice...
    vector<int> ret(ref_to_hashes.size(), 0);
    #pragma omp parallel for
    for (int i = 0; i < ref_to_hashes.size(); i++){
         vector<hash_t> matches = hash_intersection(read_hashes, ref_to_hashes[i].second);
         ret[i] = matches.size();
    }
    return ret;
    //return std::make_tuple(sample, shared_intersection, total_union);   
};

inline tuple<string, int, int> classify_and_count(vector<hash_t>& read_mins, vector<string> ref_keys, vector<vector<hash_t> > ref_mins){
    int max_shared = 0;
    string sample = "";
    int shared_inter = 0;
    int total_union = 0;
    
    for (int i = 0; i < ref_mins.size(); i++){
        vector<hash_t> matches = hash_intersection(read_mins, ref_mins[i]);
        if (matches.size() > max_shared){
            max_shared = matches.size();
            sample = ref_keys[i];
            shared_inter = matches.size();
            total_union = read_mins.size() < ref_mins[i].size() ? read_mins.size() : ref_mins[i].size();
        }
    }
    return std::make_tuple(sample, shared_inter, total_union);
};

inline tuple<string, int, int> classify_and_count(vector<string> ref_keys, vector<hash_t*> ref_mins, hash_t* read_mins,
                                                    int* ref_starts, int read_start,
                                                    int* ref_lens, int read_len,
                                                    int sketch_size){
    int max_shared = 0;
    string sample = "";
    int shared_inter = 0;
    int total_union = 0;
  
    /*
     *  std::tuple<hash_t*, int> hash_intersection(hash_t* alpha, int alpha_len, int alpha_start,
     *                                                              hash_t* beta, int beta_len, int beta_start,
     *                                                                                                                          int sketch_size){
     * */
    for (int i = 0; i < ref_keys.size(); i++){
        std::tuple<hash_t*, int> inter = hash_intersection(ref_mins[i], ref_starts[i], ref_lens[i], read_mins, read_start, read_len, sketch_size);
        int shared = std::get<1>(inter);
        if (shared > max_shared){
            max_shared = shared;
            sample = ref_keys[i];
            shared_inter = shared;
            total_union = read_len < ref_lens[i] ? read_len : ref_lens[i];
        }
        delete [] std::get<0>(inter);

    }
    return std::make_tuple(sample, shared_inter, total_union);
};



inline tuple<string, int, int, bool> classify_and_count_diff_filter(vector<string> ref_keys, vector<hash_t*> ref_mins, hash_t* read_mins,
                                                    int* ref_starts, int read_start,
                                                    int* ref_lens, int read_len,
                                                    int sketch_size, int min_diff){
    int max_shared = 0;
    int prev_best = 0;
    string sample = "";
    int shared_inter = 0;
    int total_union = 0;
  
    /*
     *  std::tuple<hash_t*, int> hash_intersection(hash_t* alpha, int alpha_len, int alpha_start,
     *                                                              hash_t* beta, int beta_len, int beta_start,
     *                                                                                                                          int sketch_size){
     * */
    for (int i = 0; i < ref_keys.size(); i++){
        std::tuple<hash_t*, int> inter = hash_intersection(ref_mins[i], ref_starts[i], ref_lens[i], read_mins, read_start, read_len, sketch_size);
        int shared = std::get<1>(inter);
        if (shared > max_shared){
            prev_best = max_shared;
            max_shared = shared;
            sample = ref_keys[i];
            shared_inter = shared;
            total_union = read_len < ref_lens[i] ? read_len : ref_lens[i];
        }
        delete [] std::get<0>(inter);

    }
    return std::make_tuple(sample, shared_inter, total_union, (max_shared - prev_best > min_diff));
};
inline tuple<string, int, int> classify_and_count_par(vector<string>& ref_keys, vector<hash_t*>& ref_mins, hash_t* read_mins,
                                                    int* ref_starts, int& read_start,
                                                    int* ref_lens, int& read_len,
                                                    int& sketch_size){
    int max_shared = 0;
    string sample = "";
    int shared_inter = 0;
    int total_union = 0;
  
    for (int i = 0; i < ref_keys.size(); i++){
        std::tuple<hash_t*, int> inter = hash_intersection(ref_mins[i], ref_starts[i], ref_lens[i], read_mins, read_start, read_len, sketch_size);
        int shared = std::get<1>(inter);
        if (shared > max_shared){
            sample = ref_keys[i];
            max_shared = shared;
            shared_inter = shared;
            total_union = read_len < ref_lens[i] ? read_len : ref_lens[i];
        }
        delete [] std::get<0>(inter);

    }
    return std::make_tuple(sample, shared_inter, total_union);
};




inline tuple<string, int, int> classify_and_count(vector<hash_t>& read_hashes, map<string, vector<hash_t> >& ref_to_hashes){
     //Parallel compare through the map would be really nice...
    int max_shared = 0;
    string sample = "";
    int shared_intersection = 0;
    int total_union = 0;
    map<string, vector<hash_t> >::iterator iter;
    for (iter = ref_to_hashes.begin(); iter != ref_to_hashes.end(); iter++){
         vector<hash_t> matches = hash_intersection(read_hashes, iter->second);
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


inline tuple<string, int, int> p_classify_and_count(vector<hash_t>& read_hashes, map<string, vector<hash_t> >& ref_to_hashes){
     //Parallel compare through the map would be really nice...
    int max_shared = 0;
    string sample = "";
    int shared_intersection = 0;
    int total_union = 0;
    vector<pair<string, vector<hash_t> > > ref_pairs(ref_to_hashes.begin(), ref_to_hashes.end());
    #pragma omp parallel for
    for (int i = 0; i < ref_pairs.size(); i++){
         vector<hash_t> matches = hash_intersection(read_hashes, ref_pairs[i].second);
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

inline vector<int> all_hash_compare(vector<hash_t>& hashes, vector<pair<string, vector<hash_t> > >& ref_hashes){
    vector<int> ret(ref_hashes.size(), 0);
    #pragma omp for
    for (int i = 0; i < ret.size(); i++){
        ret[i] = hash_intersection(hashes, ref_hashes[i].second).size();
    }

    return ret;
    
};

inline string classify(vector<hash_t>& read_hashes, map<string, vector<hash_t> >& ref_to_hashes){

    //Parallel compare through the map would be really nice...
    int max_shared = 0;
    string ret = "";
    map<string, vector<hash_t> >::iterator iter;
    for (iter = ref_to_hashes.begin(); iter != ref_to_hashes.end(); iter++){
         vector<hash_t> matches = hash_intersection(read_hashes, iter->second);
         if (matches.size() > max_shared){
            ret = iter->first;
            max_shared = matches.size();
         }
    }
    return ret;
};

#endif
