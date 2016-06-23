#include "mkmh.hpp"
#include <string>
#include <iostream>

using namespace std;
using namespace mkmh;

bool testify(int t_num, string test, string obs, string exp){
    if (obs == exp){
        cout << "PASS: " << t_num << " " << test << endl;
        return true;
    }
    else{
        cout << "FAIL: " << t_num << " " << test << endl;
        return true;
    }
}

bool testify(int t_num, string test, bool x){
    if (x){
        testify(t_num, test, "", "");
        return true;
    }
    else{
        testify(t_num, test, "1", "2");
        return false;
    }   
}

int main(){
    int t_num = 0;

    string seq = "ATGCATGCATGCATGCATGC";
    vector<string> kmers = kmerize(seq, 5);
    bool x = true;
    for (auto e : kmers){
        if (e.length() != 5){
            x = false;
        }
    }

    /* 
     *
     * Test kmerize, shingle, and minhash64
     * 
     * */
    testify(t_num++, "The kmers produced by kmerize are the right size", x);

    vector<string> shingles = shingle(seq, 5);
    x = true;
    for (auto e : shingles){
        if (e.length() != 5){
            x = false;
        }
    }
    testify(t_num++, "The shingles produced by kmerize are the right size", x);

    vector<int64_t> ret = minhash_64(seq, 5, 5);
    testify(t_num++, "minhash_64 produces the right number of hashes", ret.size() == 5);

    vector<string> k_set = kmer_set(kmers);
    testify(t_num++, "Kmer set removes duplicate kmers", k_set.size() < kmers.size());

    /* Test top_minhash_64 and bottom_minhash_64 */
    x = true;
    vector<int64_t> tops = top_minhash_64(seq, 5, 5);
    vector<int64_t> bottoms = bottom_minhash_64(seq, 5, 5);
    for (auto e : tops){
        for (auto f : bottoms){
            if (e < f){
                x = false;
            }
        }
    }
    testify(t_num++, "top_minhash64 produces the bigger values than bottom_minhash_64", x);


    string seq2 = "ACTGaaatttt";
    vector<int> ks;
    ks.push_back(4);
    ks.push_back(4);
    kmers = kmerize(seq2, 4);
    vector<string> m_kmers = multi_kmerize(seq2, ks);

    /* Test multi_kmerize */
    testify(t_num++, "multi_kmerize produces twice as many kmers with two kmer sizes", kmers.size() * 2 == m_kmers.size());

    /* Test hash union / intersection */
    vector<int64_t> t1 = {1, 2, 3, 4, 5, 6};
    vector<int64_t> t2 = {4, 5, 6, 7, 8, 9};
    testify(t_num++, "Hash intersection of two sets is the expected size.", hash_intersection(t1, t2).size() == 3);

    t1 = {1, 2, 3, 4, 4, 4, 5, 6};
    t2 = {4, 4, 5, 6, 7, 8, 9};
    testify(t_num++, "Hash intersection counts duplicate values", hash_intersection(t1, t2).size() == 4);

    testify(t_num++, "Hash union yields the correct size with duplicate values", hash_union(t1, t1).size() ==   2 * t1.size());
    //testify(t_num++ "Union / Intersection produces expected value", hash_intersection(t1, t1).size() / hash_intersection(
    testify(t_num++, "Hash union yields the correct size when duplicates are removed", v_set(hash_union(t1, t1)).size() == 6);

    testify(t_num++, "Hash set union produces the expected number of values", hash_set_union(t1, t2).size() == 9);
    testify(t_num++, "Hash set intersection produces the expected number of values", hash_set_intersection(t1, t2).size() == 3);


    return 0;
}
