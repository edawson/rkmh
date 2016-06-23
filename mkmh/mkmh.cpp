#include "mkmh.hpp"

namespace mkmh{

    using namespace std;
    
    string reverse_complement(string seq){
        stringstream ret;

        for (int i = 0; i < seq.length(); i++){
            char c = seq[i];
            switch (c){
                case 'A':
                    ret << "T";
                    break;
                case 'a':
                    ret << "c";
                    break;
                case 'T':
                    ret << "A";
                    break;
                case 't':
                    ret << "a";
                    break;
                case 'C':
                    ret << "G";
                    break;
                case 'c':
                    ret << "g";
                    break;
                case 'G':
                    ret << "C";
                    break;
                case 'g':
                    ret << "c";
                    break;
                /* Handle X, N, Y, all that stuff. */
                default:
                    ret << c;
                    break;
            }
        }
        return ret.str();

    }

    string reverse(string seq){
        string copy = string(seq);
        std::reverse(copy.begin(), copy.end());
        return copy;
    }

    vector<string> kmer_set(vector<string> kmers){
        set<string> uniqs = set<string> (kmers.begin(), kmers.end());
        vector<string> ret = vector<string> (uniqs.begin(), uniqs.end());
        return ret;
    }

    /* Returns the forward and reverse-reverse complement kmers of a sequence */
    vector<string> kmerize(string seq, int k){
        int i = 0;
        vector<string> ret;
        for (i = 0; i + k < seq.length(); i++){
            string s = seq.substr(i, k);
            ret.push_back(s);
            ret.push_back(reverse(reverse_complement(s)));
        }
        return ret;
    }

    vector<string> multi_kmerize(string seq, vector<int> kSizes){
        int i = 0;
        vector<string> ret;
        for (auto k : kSizes){
            vector<string> kmers = kmerize(seq, k);
            ret.reserve(ret.size() + kmers.size());
            ret.insert(ret.end(), kmers.begin(), kmers.end());
            //for (i = 0; i + k < seq.length(); i++){
            //    ret.push_back(seq.substr(i, i+k));
            //    ret.push_back(reverse(reverse_complement(seq.substr(i, i+k))));
            //}
        }
        return ret;
    }

    /* Returns the forward shingles size k of a sequence */
    vector<string> shingle(string seq, int k){
        int i = 0;
        vector<string> ret;
        for (i = 0; i < seq.length() - k; i++){
            ret.push_back(seq.substr(i, k));
        }
        return ret;
    }

    vector<string> multi_shingle(string seq, vector<int> kSizes){
        int i = 0;
        vector<string> ret;
        for (auto k : kSizes){
            for (i = 0; i + k < seq.length(); i++){
                ret.push_back(seq.substr(i, k));
            }
        }
        return ret;
    }

    // vector<int64_t> preserve_kmer_mh64(string seq, vector<int> kSizes, int hashSize);

    vector<int64_t> minhash_64(string seq, int k, int hashSize, bool useBottom){
        vector<int64_t> ret;
        vector<string> kmers = kmerize(seq, k);
        ret.reserve(kmers.size());
        uint32_t seed = 101;

        vector<string>::iterator it;
        for (it = kmers.begin(); it != kmers.end(); it++){
            uint32_t khash[4];
            MurmurHash3_x64_128(&(*it), (*it).length(), seed, khash);
            //MurmurHash3_x64_128(argv[1], strlen(argv[1]), seed, hash);
            int64_t r_hash = int64_t(khash[2]) << 32 | int64_t(khash[1]);
            ret.push_back(r_hash);
        }

        std::sort(ret.begin(), ret.end());

        return useBottom ?
            vector<int64_t> (ret.begin(), ret.begin() + hashSize) :
            vector<int64_t> (ret.end() - hashSize ,ret.end());
    }

    vector<int64_t> top_minhash_64(string seq, int k, int hashSize){
        return minhash_64(seq, k, hashSize, false);
    }

    vector<int64_t> bottom_minhash_64(string seq, int k, int hashSize){
        return minhash_64(seq, k, hashSize, true);
    }

    vector<int64_t> hash_intersection(vector<int64_t> alpha, vector<int64_t> beta){
        vector<int64_t> ret;
        int i = 0;
        int j = 0;
        while (i < alpha.size() && j < beta.size()){
            if (alpha[i] == beta[j]){
                ret.push_back(alpha[i]);
                i++;
                j++;
            }
            else if (alpha[i] > beta[j]){
                j++;
            }
            else{
                i++;
            }
        }

        return ret;
    }

    vector<int64_t> hash_union(vector<int64_t> alpha, vector<int64_t> beta){
        vector<int64_t> ret;


        ret.reserve(alpha.size() + beta.size());
        ret = vector<int64_t> (alpha.begin(), alpha.end());
        ret.insert(ret.end(), beta.begin(), beta.end());
        return ret;
    }

    vector<int64_t> hash_set_intersection(vector<int64_t> alpha, vector<int64_t> beta){
        return hash_intersection(v_set(alpha), v_set(beta));
        
    }

    vector<int64_t> hash_set_union(vector<int64_t> alpha, vector<int64_t> beta){
        return v_set(hash_union(alpha, beta));
    }

}
