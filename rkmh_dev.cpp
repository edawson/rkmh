#include <iostream>
#include <algorithm>
#include <vector>
#include<set>
#include <cstdint>
#include <string>
#include <list>
#include <sstream>
#include <zlib.h>
#include <omp.h>
#include <getopt.h>
#include <map>
#include <unordered_map>
#include "mkmh.hpp"
#include "kseq.hpp"
#include "equiv.hpp"

KSEQ_INIT(gzFile, gzread)

    using namespace std;
    using namespace mkmh;

    /**
     * Proposed CLI:
     * ./rkmh classify
     * sketch size: s
     * kmer: k
     * readfile: f
     * reference file: r
     * thread: t
     * min_kmer_occ: M
     * minmatches: N
     * min diff: D
     * write to outfile: o
     * min informative: I
     * pre-calculated hashes: p
     * call differing variants: c
     *
     * ./rkmh hash
     * sketch size: -s
     * min_kmer_occ: -M
     * min_samples: -I
     * input: -f
     *
     * ./rkmh call
     * same as classify
     * kmer size: k
     * readfile: f
     * ref file: r
     * min_inform: I
     * min diff: D
     * kmer_occ: M
     * outfile: o
     * precalc: p
     * minimum depth: -d
     */

    void print_help(char** argv){
        cerr << "Usage: " << argv[0] << " [options]" << endl
            << "Options:" << endl
            << "--reference/-r   <REF>" << endl
            << "--fasta/-f   <FASTAFILE>" << endl
            << "--kmer/-k    <KMERSIZE>" << endl
            << "--sketch-size/-s <SKETCHSIZE>" << endl
            << "--threads/-t <THREADS>" << endl
            << "--min-kmer-occurence/-M <MINOCCURENCE>" << endl
            << "--min-matches/-N <MINMATCHES>" << endl
            << "--min-diff/-D    <MINDIFFERENCE>" << endl
            << "--min-informative/-I <MAXSAMPLES> only use kmers present in fewer than MAXSAMPLES" << endl
            << endl;
    }

void parse_fastas(vector<char*>& files,
        unordered_map<string, char*>& ret_to_seq,
        unordered_map<string, int>& ret_to_len){

    /*
    gzFile fp;
    kseq_t *seq;
    int l;
    */

    for (auto f : files){
        gzFile fp;
        kseq_t *seq;
        int l;
        fp = gzopen(f, "r");
        seq = kseq_init(fp);
        // Read in reads, cluster, spit it back out
        while ((l = kseq_read(seq)) >= 0) {
            to_upper(seq->seq.s, seq->seq.l);

            char * x = new char[seq->seq.l];
            memcpy(x, seq->seq.s, seq->seq.l);
            ret_to_seq[string(seq->name.s)] = x; 
            ret_to_len[seq->name.s] = seq->seq.l;
        } 
        gzclose(fp);
    }
}

void parse_fastas(vector<char*>& files,
                vector<string>& seq_keys,
                vector<char*>& seq_seqs,
                vector<int>& seq_lens){

    kseq_t *seq;
     for (int i = 0; i < files.size(); i++){
        char* f = files[i];
        gzFile fp;
        int l;
        fp = gzopen(f, "r");
        seq = kseq_init(fp);
        // Read in reads, cluster, spit it back out
        while ((l = kseq_read(seq)) >= 0) {
            to_upper(seq->seq.s, seq->seq.l);

            char * x = new char[seq->seq.l];
            memcpy(x, seq->seq.s, seq->seq.l);
            seq_keys.push_back(seq->name.s);
            seq_seqs.push_back(x);
            seq_lens.push_back(seq->seq.l);
        } 
        gzclose(fp);
    }   
    kseq_destroy(seq);
}

void hash_sequences(vector<string>& keys,
    unordered_map<string, char*>& name_to_seq,
    unordered_map<string, int>& name_to_length,
    vector<int>& kmer,
    unordered_map<string, hash_t*>& ret_to_hashes,
    unordered_map<string, int>& ret_to_hash_num){

    #pragma omp for
    for (int i = 0; i < keys.size(); i++){
        tuple<hash_t*, int> hashes_and_num =  allhash_unsorted_64_fast(name_to_seq[keys[i]], name_to_length[keys[i]], kmer);
        ret_to_hashes[keys[i]] = std::get<0>(hashes_and_num);
        ret_to_hash_num[keys[i]] = std::get<1>(hashes_and_num);
    }

}


/**
 * Requires that all vectors be of the same length
 */
void hash_sequences(vector<string>& keys,
    vector<char*>& seqs,
    vector<int>& lengths,
    vector<hash_t*>& hashes,
    vector<int>& hash_lengths,
    vector<int>& kmer,
    unordered_map<hash_t, int>& read_hash_to_depth,
    unordered_map<hash_t, int>& ref_to_sample_depth,
    bool doReadDepth,
    bool doReferenceDepth){

    #pragma omp parallel for
    for (int i = 0; i < keys.size(); i++){
        // Hash sequence
        tuple<hash_t*, int> hashes_and_num =  allhash_unsorted_64_fast(seqs[i], lengths[i], kmer);
        hashes[i] = std::get<0>(hashes_and_num);
        hash_lengths[i] = std::get<1>(hashes_and_num);
        if (doReadDepth){
            #pragma omp critical
            {
            for (int j = 0; j < hash_lengths[i]; j++){
                #pragma omp atomic update
                read_hash_to_depth[hashes[i][j]] ++;
            }
            }
        }
        else if (doReferenceDepth){
           // create the set of hashes in the sample
           set<hash_t> sample_set (hashes[i], hashes[i] + hash_lengths[i]);
            #pragma omp critical
           {
           for (auto x : sample_set){
                //#pragma omp atomic update 
                ref_to_sample_depth[x] ++;
           }
           }
        }
        
    }

}




/**
 *
 *     * ./rkmh call
 * same as classify
 * kmer size: k
 * readfile: f
 * ref file: r
 * min_inform: I
 * min diff: D
 * kmer_occ: M
 * outfile: o
 * precalc: p
 * minimum depth: -d
 *
 */
int main_call(int argc, char** argv){
    vector<char*> ref_files;
    vector<char*> read_files;


    unordered_map<hash_t, int> read_hash_to_depth;
    read_hash_to_depth.reserve(10000);
    unordered_map<hash_t, int> ref_hash_to_num_samples;
    ref_hash_to_num_samples.reserve(10000);

    vector<int> kmer;

    int sketch_size = 1000;
    int threads = 1;
    int min_kmer_occ = 1;
    int min_matches = -1;
    int min_diff = 0;
    int max_samples = 1000000;

    int c;
    int optind = 2;

    if (argc <= 2){
        print_help(argv);
        exit(1);
    }

    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"kmer", no_argument, 0, 'k'},
            {"fasta", required_argument, 0, 'f'},
            {"reference", required_argument, 0, 'r'},
            {"sketch", required_argument, 0, 's'},
            {"threads", required_argument, 0, 't'},
            {"min-kmer-occurence", required_argument, 0, 'M'},
            {"min-matches", required_argument, 0, 'N'},
            {"min-diff", required_argument, 0, 'D'},
            {"max-samples", required_argument, 0, 'I'},
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "hk:f:r:s:t:M:N:D:I:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case 't':
                threads = atoi(optarg);
                break;
            case 'r':
                ref_files.push_back(optarg);
                break;
            case 'f':
                read_files.push_back(optarg);
                break;
            case 'k':
                kmer.push_back(atoi(optarg));
                break;
            case '?':
            case 'h':
                print_help(argv);
                exit(1);
                break;
            case 's':
                sketch_size = atoi(optarg);
                break;
            case 'M':
                min_kmer_occ = atoi(optarg);
                break;
            case 'N':
                min_matches = atoi(optarg);
                break;
            case 'D':
                min_diff = atoi(optarg);
                break;
            case 'I':
                max_samples = atoi(optarg);
                break;
            default:
                print_help(argv);
                abort();

        }
    }

    omp_set_num_threads(threads);

    vector<string> ref_keys;
    ref_keys.reserve(500);
    vector<char*> ref_seqs;
    ref_seqs.reserve(500);
    vector<int> ref_lens;
    ref_lens.reserve(500);
    vector<string> read_keys;
    read_keys.reserve(2000);
    vector<char*> read_seqs;
    read_seqs.reserve(2000);
    vector<int> read_lens;
    read_lens.reserve(2000);

    #pragma omp master
    cerr << "Parsing sequences...";
    parse_fastas(ref_files, ref_keys, ref_seqs, ref_lens);
    parse_fastas(read_files, read_keys, read_seqs, read_lens);
    
    vector<hash_t*> ref_hashes(ref_keys.size());
    vector<int> ref_hash_nums(ref_keys.size());

    vector<hash_t*> read_hashes(read_keys.size());
    vector<int> read_hash_nums(read_keys.size());



//#pragma omp parallel
{
    
    #pragma omp master
    cerr << " Done." << endl <<
        ref_keys.size() << " references and " << read_keys.size() << " reads parsed." << endl;
    /*
     *void hash_sequences(vector<string>& keys,
    vector<char*>& seqs,
    vector<int>& lengths,
    vector<hash_t*>& hashes,
    vector<int>& hash_lengths,
    unordered_map<hash_t, int>& read_hash_to_depth,
    unordered_map<hash_t, int>& ref_to_sample_depth,
    bool doReadDepth,
    bool doReferenceDepth){
     */
   
#pragma omp master
    cerr << "Hashing references... ";
    hash_sequences(ref_keys, ref_seqs, ref_lens,
                  ref_hashes, ref_hash_nums, kmer,
                  read_hash_to_depth,
                  ref_hash_to_num_samples,
                  (min_kmer_occ > 0),
                  (max_samples < 10000));
#pragma omp master
    cerr << " Done." << endl;

#pragma omp master
    cerr << "Hashing reads... ";
    hash_sequences(read_keys, read_seqs, read_lens,
                    read_hashes, read_hash_nums, kmer,
                    read_hash_to_depth,
                    ref_hash_to_num_samples,
                    (min_kmer_occ > 0),
                    (max_samples < 10000));
    #pragma omp master
    cerr << " Done." << endl;
    vector<vector<hash_t> > ref_mins(ref_keys.size(), vector<hash_t>(1));

    /**
     *
     * vector<char> ATGC = {'A', 'T', 'G', 'C'};
     * while (i % 4 != y)
     * auto permute = [](char* x, int k){
     *
     *  snps: switch-case
     *  for (int i = 0; i < k; i++){
     *      char val = x[i];
     *      switch (val){
     *          case 'A':
     *          case 'a':
     *              
     *
     *      }
     *  }
     *
     *  indels: for loops
     * }
     */
    list<int> d_window;
    auto avg = [](list<int> n_list){
        int ret = 0;
        //for (int x = 0; x < n_list.size(); x++){
        //    ret += ;
        //}
        return (double) ret / (double) n_list.size();
    };
    #pragma omp for
    for (int i = 0; i < ref_keys.size(); i++){
        for (int j = 0; j < ref_hash_nums[i]; j++){
            int depth = read_hash_to_depth[ref_hashes[i][j]];
            d_window.push_back(depth);
            if (d_window.size() > 5){
                d_window.pop_front();
            }
            if (depth < .5 * avg(d_window)){

            }

        }
        
    }

   }

    for (auto x : read_hashes){
        delete [] x;;
    }
    for (auto x : read_seqs){
        delete [] x;
    }
    for (auto y : ref_hashes){
        delete [] y;
    }
    for (auto y : ref_seqs){
        delete [] y;
    }

    return 0;
}
    // Parse CLI
    //
    // Generate hashes and filter
    //parse_fastas(ref_files, ref_keys, ref_seqs, ref_lens);
    //parse_fasta(ref_files, ref_keys, ref_seqs, ref_lens);
    //
    //hash_sequences(ref_keys,
    //
    //hash_sequences(ref_keys,
    //
    //
    // Classify reads
    //
    // Call variants
    //      calculate avg depth
    //      check depth
    //
    //      permute ref sequence
    //      check new depth
    //      emit calls (if any)



/**
 *
 */
int main_hash(int argc, char** argv){
    vector<char*> ref_files;
    vector<char*> read_files;

    unordered_map<hash_t, int> read_hash_to_depth;
    read_hash_to_depth.reserve(10000);
    unordered_map<hash_t, int> ref_hash_to_num_samples;
    ref_hash_to_num_samples.reserve(10000);

    vector<int> kmer;

    int sketch_size = 1000;
    int threads = 1;
    int min_kmer_occ = 0;
    int min_matches = -1;
    int min_diff = 0;
    int max_samples = 1000000;

    int c;
    int optind = 2;

    if (argc <= 2){
        print_help(argv);
        exit(1);
    }

    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"kmer", no_argument, 0, 'k'},
            {"fasta", required_argument, 0, 'f'},
            {"reference", required_argument, 0, 'r'},
            {"sketch", required_argument, 0, 's'},
            {"threads", required_argument, 0, 't'},
            {"min-kmer-occurence", required_argument, 0, 'M'},
            {"min-matches", required_argument, 0, 'N'},
            {"min-diff", required_argument, 0, 'D'},
            {"max-samples", required_argument, 0, 'I'},
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "hk:f:r:s:t:M:N:D:I:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case 't':
                threads = atoi(optarg);
                break;
            case 'r':
                ref_files.push_back(optarg);
                break;
            case 'f':
                read_files.push_back(optarg);
                break;
            case 'k':
                kmer.push_back(atoi(optarg));
                break;
            case '?':
            case 'h':
                print_help(argv);
                exit(1);
                break;
            case 's':
                sketch_size = atoi(optarg);
                break;
            case 'M':
                min_kmer_occ = atoi(optarg);
                break;
            case 'N':
                min_matches = atoi(optarg);
                break;
            case 'D':
                min_diff = atoi(optarg);
                break;
            case 'I':
                max_samples = atoi(optarg);
                break;
            default:
                print_help(argv);
                abort();

        }
    }

    omp_set_num_threads(threads);

    vector<string> ref_keys;
    ref_keys.reserve(500);
    vector<char*> ref_seqs;
    ref_seqs.reserve(500);
    vector<int> ref_lens;
    ref_lens.reserve(500);
    vector<string> read_keys;
    read_keys.reserve(2000);
    vector<char*> read_seqs;
    read_seqs.reserve(2000);
    vector<int> read_lens;
    read_lens.reserve(2000);

    #pragma omp master
    cerr << "Parsing sequences...";
    parse_fastas(ref_files, ref_keys, ref_seqs, ref_lens);
    parse_fastas(read_files, read_keys, read_seqs, read_lens);
    
    vector<hash_t*> ref_hashes(ref_keys.size());
    vector<int> ref_hash_nums(ref_keys.size());

    vector<hash_t*> read_hashes(read_keys.size());
    vector<int> read_hash_nums(read_keys.size());



#pragma omp parallel
{
    
    #pragma omp master
    cerr << " Done." << endl <<
        ref_keys.size() << " references and " << read_keys.size() << " reads parsed." << endl;

    #pragma omp master
    cerr << "Hashing references... ";
    hash_sequences(ref_keys, ref_seqs, ref_lens,
                  ref_hashes, ref_hash_nums, kmer,
                  read_hash_to_depth,
                  ref_hash_to_num_samples,
                  (min_kmer_occ > 0),
                  (max_samples < 10000));
#pragma omp master
    cerr << " Done." << endl;

#pragma omp master
    cerr << "Hashing reads... ";
    hash_sequences(read_keys, read_seqs, read_lens,
                    read_hashes, read_hash_nums, kmer,
                    read_hash_to_depth,
                    ref_hash_to_num_samples,
                    (min_kmer_occ > 0),
                    (max_samples < 10000));
    #pragma omp master
    cerr << " Done." << endl;
    vector<vector<hash_t> > ref_mins(ref_keys.size(), vector<hash_t>(1));

    // min refs

    //min reads

}


    
       return 0;
}


/**
 *
 * ./rkmh classify
     * sketch size: s
     * kmer: k
     * readfile: f
     * reference file: r
     * thread: t
     * min_kmer_occ: M
     * minmatches: N
     * min diff: D
     * write to outfile: o
     * min informative: I
     * pre-calculated hashes: p
     * call differing variants: c

 */

int main_classify(int argc, char** argv){

    vector<char*> ref_files;
    vector<char*> read_files;


    unordered_map<hash_t, int> read_hash_to_depth;
    read_hash_to_depth.reserve(10000);
    unordered_map<hash_t, int> ref_hash_to_num_samples;
    ref_hash_to_num_samples.reserve(10000);

    vector<int> kmer;

    int sketch_size = -1;
    int threads = 1;
    int min_kmer_occ = 0;
    int min_matches = -1;
    int min_diff = 0;
    int max_samples = 1000000;

    int c;
    int optind = 2;

    if (argc <= 2){
        print_help(argv);
        exit(1);
    }

    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"kmer", no_argument, 0, 'k'},
            {"fasta", required_argument, 0, 'f'},
            {"reference", required_argument, 0, 'r'},
            {"sketch", required_argument, 0, 's'},
            {"threads", required_argument, 0, 't'},
            {"min-kmer-occurence", required_argument, 0, 'M'},
            {"min-matches", required_argument, 0, 'N'},
            {"min-diff", required_argument, 0, 'D'},
            {"max-samples", required_argument, 0, 'I'},
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "hk:f:r:s:t:M:N:D:I:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case 't':
                threads = atoi(optarg);
                break;
            case 'r':
                ref_files.push_back(optarg);
                break;
            case 'f':
                read_files.push_back(optarg);
                break;
            case 'k':
                kmer.push_back(atoi(optarg));
                break;
            case '?':
            case 'h':
                print_help(argv);
                exit(1);
                break;
            case 's':
                sketch_size = atoi(optarg);
                break;
            case 'M':
                min_kmer_occ = atoi(optarg);
                break;
            case 'N':
                min_matches = atoi(optarg);
                break;
            case 'D':
                min_diff = atoi(optarg);
                break;
            case 'I':
                max_samples = atoi(optarg);
                break;
            default:
                print_help(argv);
                abort();

        }
    }

    omp_set_num_threads(threads);

    vector<string> ref_keys;
    ref_keys.reserve(500);
    vector<char*> ref_seqs;
    ref_seqs.reserve(500);
    vector<int> ref_lens;
    ref_lens.reserve(500);
    vector<string> read_keys;
    read_keys.reserve(2000);
    vector<char*> read_seqs;
    read_seqs.reserve(2000);
    vector<int> read_lens;
    read_lens.reserve(2000);

    #pragma omp master
    cerr << "Parsing sequences...";
    parse_fastas(ref_files, ref_keys, ref_seqs, ref_lens);
    parse_fastas(read_files, read_keys, read_seqs, read_lens);
    
    vector<hash_t*> ref_hashes(ref_keys.size());
    vector<int> ref_hash_nums(ref_keys.size());

    vector<hash_t*> read_hashes(read_keys.size());
    vector<int> read_hash_nums(read_keys.size());



//#pragma omp parallel
{
    
    #pragma omp master
    cerr << " Done." << endl <<
        ref_keys.size() << " references and " << read_keys.size() << " reads parsed." << endl;
    /*
     *void hash_sequences(vector<string>& keys,
    vector<char*>& seqs,
    vector<int>& lengths,
    vector<hash_t*>& hashes,
    vector<int>& hash_lengths,
    unordered_map<hash_t, int>& read_hash_to_depth,
    unordered_map<hash_t, int>& ref_to_sample_depth,
    bool doReadDepth,
    bool doReferenceDepth){
     */
   
#pragma omp master
    cerr << "Hashing references... ";
    hash_sequences(ref_keys, ref_seqs, ref_lens,
                  ref_hashes, ref_hash_nums, kmer,
                  read_hash_to_depth,
                  ref_hash_to_num_samples,
                  false,
                  (max_samples < 10000));
#pragma omp master
    cerr << " Done." << endl;

#pragma omp master
    cerr << "Hashing reads... ";
    hash_sequences(read_keys, read_seqs, read_lens,
                    read_hashes, read_hash_nums, kmer,
                    read_hash_to_depth,
                    ref_hash_to_num_samples,
                    (min_kmer_occ > 0),
                    false);
    #pragma omp master
    cerr << " Done." << endl;
    vector<vector<hash_t> > ref_mins(ref_keys.size(), vector<hash_t>(1));
    #pragma omp parallel for
    for (int i = 0; i < ref_keys.size(); i++){
        vector<hash_t> x;
        std::sort(ref_hashes[i], ref_hashes[i] + ref_hash_nums[i]);
        if (max_samples < 10000){
            for (int j = 0; j < ref_hash_nums[i]; j++){
                if (x.size() >= sketch_size){
                    ref_mins[i] = x;
                    break;
                }
                if (ref_hashes[i][j] == 0){
                    continue;
                }
                else if (ref_hash_to_num_samples[ref_hashes[i][j]] > max_samples){
                    continue;
                }
                else{
                    x.push_back(ref_hashes[i][j]);
                }
            }
            
            //map<string, vector<hash_t> > informs = only_informative_kmers(ref_to_hashes, ref_to_num_hashes, max_samples);
            //for (int i = 0; i < ref_seq.size(); i++){
            //    vector<hash_t> mins = minhashes(&(informs[ref_seq[i].first]).front(), informs[ref_seq[i].first].size(), sketch_size);
            //    ref_to_mins[ref_seq[i].first] = mins;
            //}
            // set ref_to_hashes to a new value
            // set ref_to_num_hashes to a new value
            // need a only_informative_kmers(hash_t*, int hash_num, int max_sample) function
        }
        else{
            x = minhashes(ref_hashes[i], ref_hash_nums[i], sketch_size);
            ref_mins[i] = x;
        }



    }

    #pragma omp parallel for
    for (int i = 0; i < read_keys.size(); i++){
        hash_t* hh = read_hashes[i];
        int hh_l = read_hash_nums[i];
        
        //generate mins

        // compare / classify / print to cout

            stringstream outre;

            sort(hh, hh + hh_l);
            vector<hash_t> mins;
            if (min_kmer_occ > 0){
                    int ret_ind;
                    while (hh[ret_ind] == 0){
                        ret_ind++;
                    }
                    int hashmax = sketch_size + ret_ind < hh_l ? sketch_size + ret_ind : hh_l - 1;

                    for (int j = 0; j < hh_l; j++){
                        if (mins.size() == hashmax ){
                            break;
                        }

                        if (read_hash_to_depth[hh[j]] >= min_kmer_occ && hh[j] != 0){
                            mins.push_back(hh[j]);
                        }
                    }
                }
            else{

                /**sort(hashes, hashes + hash_len);


                  int ret_ind;
                  while (hashes[ret_ind] == 0){
                  ret_ind++;
                  }
                  int hashmax = sketch_size + ret_ind < hash_len ? sketch_size + ret_ind : hash_len - 1;

                  mins = vector<hash_t>(hashes + ret_ind, hashes + hashmax);
                 **/
                mins = minhashes(hh, hh_l, sketch_size);
            }
            tuple<string, int, int> result;
            result = classify_and_count(mins, ref_keys, ref_mins);


            bool depth_filter = mins.size() == 0; 
            bool match_filter = std::get<1>(result) < min_matches;

            outre  << "Sample: " << read_keys[i] << "\t" << "Result: " << 
                std::get<0>(result) << "\t" << std::get<1>(result) << "\t" << std::get<2>(result) << "\t" <<
                (depth_filter ? "FAIL:DEPTH" : "") << "\t" << (match_filter ? "FAIL:MATCHES" : "") << endl;

            cout << outre.str();

    }
    }

    for (auto x : read_hashes){
        delete [] x;;
    }
    for (auto x : read_seqs){
        delete [] x;
    }
    for (auto y : ref_hashes){
        delete [] y;
    }
    for (auto y : ref_seqs){
        delete [] y;
    }

    return 0;
}


int main(int argc, char** argv){
    string cmd = argv[1];
    if (cmd == "classify"){
        return main_classify(argc, argv);
    }
    else if (cmd == "hash"){

    }
    else if (cmd == "call"){

    }

    vector<char*> ref_files;
    vector<char*> read_files;
    vector<int> kmer;

    int sketch_size = -1;
    int threads = 1;
    int min_kmer_occ = 0;
    int min_matches = -1;
    int min_diff = 0;
    int max_samples = 1000000;

    stringstream errtre;

    unordered_map<string, hash_t* > read_to_hashes;
    unordered_map<string, int> read_to_num_hashes;
    unordered_map<string, char*> read_to_seq;
    unordered_map<string, int> read_to_length;
    unordered_map<string, vector<hash_t> > read_to_mins;

    unordered_map<hash_t, int> read_hash_to_depth;


    map<string, hash_t*> ref_to_hashes;
    map<string, int> ref_to_num_hashes;
    unordered_map<string, char*> ref_to_seq;
    unordered_map<string, int> ref_to_length;
    map<string, vector<hash_t> > ref_to_mins;

    unordered_map<hash_t, int> ref_hash_to_depth;

    unordered_map<string, unordered_map<int, hash_t> > ref_to_pos_to_depth;

    bool do_call = false;
    bool do_hash = false;
    bool do_classify = true;

    //hash_t d_arr [INT64_MAX]; array is too large


    int c;
    if (argc < 2){
        print_help(argv);
        exit(1);
    }

    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"kmer", no_argument, 0, 'k'},
            {"reads", required_argument, 0, 'r'},
            {"fasta", required_argument, 0, 'f'},
            {"minhash", required_argument, 0, 'm'},
            {"threads", required_argument, 0, 't'},
            {"min-kmer-occurence", required_argument, 0, 'D'},
            {"min-matches", required_argument, 0, 'S'},
            {"min-diff", required_argument, 0, 'P'},
            {"min-samples", required_argument, 0, 'I'},
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "hm:k:r:f:t:D:S:P:I:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case 't':
                threads = atoi(optarg);
                break;
            case 'f':
                //ref_file = optarg;
                ref_files.push_back(optarg);
                break;
            case 'r':
                //read_file = optarg;
                read_files.push_back(optarg);
                break;
            case 'k':
                kmer.push_back(atoi(optarg));
                break;
            case '?':
            case 'h':
                print_help(argv);
                exit(1);
                break;
            case 'm':
                sketch_size = atoi(optarg);
                break;
            case 'D':
                min_kmer_occ = atoi(optarg);
                break;
            case 'S':
                min_matches = atoi(optarg);
                break;
            case 'P':
                min_diff = atoi(optarg);
                break;
            case 'I':
                max_samples = atoi(optarg);
                break;
            default:
                print_help(argv);
                abort();

        }
    }

    omp_set_num_threads(threads);
    // Read in fastas
    gzFile fp;
    kseq_t *seq;
    int l;

    char* uppered;

    if (ref_files.size() > 0){
        for (auto f : ref_files){
            fp = gzopen(f, "r");
            seq = kseq_init(fp);
            // Read in reads, cluster, spit it back out
            while ((l = kseq_read(seq)) >= 0) {
                to_upper(seq->seq.s, seq->seq.l);

                char * x = new char[seq->seq.l];
                memcpy(x, seq->seq.s, seq->seq.l);
                ref_to_seq[string(seq->name.s)] = x; 
                ref_to_length[seq->name.s] = seq->seq.l;
            } 
            gzclose(fp);
        }

        cerr << "Loaded " << ref_to_seq.size() << " reference sequences." << endl;
    }



    else{
        cerr << "Please provide a fasta file containing references." << endl;
        exit(1);
    }


    if (read_files.size() > 0){
        for (auto f : read_files){
            fp = gzopen(f, "r");
            seq = kseq_init(fp);
            while ((l = kseq_read(seq)) >= 0) {
                to_upper(seq->seq.s, seq->seq.l);

                char * x = new char[seq->seq.l];
                memcpy(x, seq->seq.s, seq->seq.l);
                read_to_seq[string(seq->name.s)] = x; 
                read_to_length[seq->name.s] = seq->seq.l;
            }
            gzclose(fp);
        }

#pragma omp master
        cerr << "Loaded " << read_to_seq.size() << " reads." << endl;
    }
    else{
        cerr << "Please provide a read file containing query sequences." << endl;
        exit(1);
    }

    if (kmer.size() < 1){
#pragma omp master
        {
            cerr << "No kmer size provided. The default of 16 will be used." << endl;
            kmer.push_back(16);
        }
    }

    vector<pair<string, const char* > > read_seq(read_to_seq.begin(), read_to_seq.end());
    vector<pair<string, const char* > > ref_seq(ref_to_seq.begin(), ref_to_seq.end());


#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < ref_seq.size(); i++){

            tuple<hash_t*, int> hashes_and_num =  allhash_unsorted_64_fast(ref_seq[i].second, ref_to_length[ref_seq[i].first], kmer);
            hash_t* hashes = std::get<0>(hashes_and_num);
            int num_hashes = std::get<1>(hashes_and_num);

            ref_to_hashes[ref_seq[i].first] =  hashes;
            ref_to_num_hashes[ref_seq[i].first] = num_hashes;

            vector<hash_t> x = minhashes(hashes, num_hashes, sketch_size);
            ref_to_mins[ref_seq[i].first] = x;
        }

#pragma omp master
        {
            errtre << "Hashed " << ref_to_hashes.size() << " references." << endl;
            cerr << errtre.str();
            errtre.str("");
        }

        if (max_samples < 1000000){
#pragma omp master
            cerr << "Filtering common kmers in references..." << endl;

            map<string, vector<hash_t> > informs = only_informative_kmers(ref_to_hashes, ref_to_num_hashes, max_samples);
            for (int i = 0; i < ref_seq.size(); i++){
                vector<hash_t> mins = minhashes(&(informs[ref_seq[i].first]).front(), informs[ref_seq[i].first].size(), sketch_size);
                ref_to_mins[ref_seq[i].first] = mins;
            }
            // set ref_to_hashes to a new value
            // set ref_to_num_hashes to a new value
            // need a only_informative_kmers(hash_t*, int hash_num, int max_sample) function
        }


#pragma omp master
        cerr << "Hashing reads..." << endl;

#pragma omp for
        for (int i = 0; i < read_seq.size(); i++){
            tuple<hash_t*, int> hashes_and_num =  allhash_unsorted_64_fast(read_seq[i].second, read_to_length[read_seq[i].first], kmer);
            hash_t* hashes = std::get<0>(hashes_and_num);
            int num_hashes = std::get<1>(hashes_and_num);

            read_to_hashes[read_seq[i].first] = hashes;
            read_to_num_hashes[read_seq[i].first] = num_hashes;

            if (min_kmer_occ > 0){

#pragma omp critical
                {
                    for (int j = 0; j < num_hashes; j++){
                        read_hash_to_depth[ hashes[j] ] ++;
                    }
                }
            }
        }

#pragma omp master
        cerr << "Classifying reads..." << endl;

#pragma omp for
        for (int i = 0; i < read_seq.size(); i++){
            stringstream outre;

            hash_t* hashes = read_to_hashes[read_seq[i].first];
            int hash_len = read_to_num_hashes[read_seq[i].first];

            vector<hash_t> mins;
            mins.reserve(sketch_size);
            if (min_kmer_occ > 0){
#pragma omp critical
                {

                    sort(hashes, hashes + hash_len);

                    int ret_ind;
                    while (hashes[ret_ind] == 0){
                        ret_ind++;
                    }
                    int hashmax = sketch_size + ret_ind < hash_len ? sketch_size + ret_ind : hash_len - 1;

                    for (int j = 0; j < hash_len; j++){
                        if (mins.size() == hashmax ){
                            break;
                        }

                        if (read_hash_to_depth[hashes[j]] >= min_kmer_occ && hashes[j] != 0){
                            mins.push_back(hashes[j]);
                        }
                    }
                }


            }


            else{

                /**sort(hashes, hashes + hash_len);


                  int ret_ind;
                  while (hashes[ret_ind] == 0){
                  ret_ind++;
                  }
                  int hashmax = sketch_size + ret_ind < hash_len ? sketch_size + ret_ind : hash_len - 1;

                  mins = vector<hash_t>(hashes + ret_ind, hashes + hashmax);
                 **/
                mins = minhashes(hashes, hash_len, sketch_size);
            }
            tuple<string, int, int> result;
            result = classify_and_count(mins, ref_to_mins);


            bool depth_filter = mins.size() == 0; 
            bool match_filter = std::get<1>(result) < min_matches;

            outre  << "Sample: " << read_seq[i].first << "\t" << "Result: " << 
                std::get<0>(result) << "\t" << std::get<1>(result) << "\t" << std::get<2>(result) << "\t" <<
                (depth_filter ? "FAIL:DEPTH" : "") << "\t" << (match_filter ? "FAIL:MATCHES" : "") << endl;

            cout << outre.str();
        }

    }
    
    /*
    if (argc == 1){

    }

    string command = argv[1];
    if (command == "classify"){

    }
    else if (command == "call"){

    }
    else if (command == "hash"){

    }
    else{

    }
    */

    for (auto x : read_to_hashes){
        delete [] x.second;
        delete [] read_to_seq[x.first];
    }
    for (auto y : ref_to_hashes){
        delete [] y.second;
        delete [] ref_to_seq[y.first];
    }

#pragma omp master
    cerr << "Done." << endl;

    kseq_destroy(seq);
    gzclose(fp);



    return 0;
}
