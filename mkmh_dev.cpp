#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdint>
#include <string>
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
     * ./rkmh hash
     *
     * ./rkmh call
     *
     */

    void print_help(char** argv){
        cerr << "Usage: " << argv[0] << " [options]" << endl
            << "Options:" << endl
            << "--reads/-r   <READFILE>" << endl
            << "--fasta/-f   <FASTAFILE>" << endl
            << "--kmer/-k    <KMERSIZE>" << endl
            << "--minhash/-m <HASHSIZE>" << endl
            << "--threads/-t <THREADS>" << endl
            << "--min-kmer-occurence/-D <MINOCCURENCE>" << endl
            << "--min-matches/-S <MINMATCHES>" << endl
            << "--min-diff/-P    <MINDIFFERENCE>" << endl
            << "--min-informative/-I <MAXSAMPLES> only use kmers present in fewer than MAXSAMPLES" << endl
            << endl;
    }


int main(int argc, char** argv){
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
    unordered_map<string, const char*> read_to_seq;
    unordered_map<string, int> read_to_length;
    map<string, vector<hash_t> > read_to_mins;

    map<hash_t, int> read_hash_to_depth;


    map<string, hash_t*> ref_to_hashes;
    map<string, int> ref_to_num_hashes;
    map<string, const char*> ref_to_seq;
    map<string, int> ref_to_length;
    map<string, vector<hash_t> > ref_to_mins;

    unordered_map<hash_t, int> ref_hash_to_depth;

    unordered_map<string, unordered_map<int, hash_t> > ref_to_pos_to_depth;

    bool do_call = false;
    unordered_map<string, unordered_map<int, string> > ref_to_pos_to_call;

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

    //map<string, std::tuple<char*, int> > name_to_seq_len;

    if (ref_files.size() > 0){
        for (auto f : ref_files){
            fp = gzopen(f, "r");
            seq = kseq_init(fp);
            // Read in reads, cluster, spit it back out
            while ((l = kseq_read(seq)) >= 0) {
                ref_to_seq[seq->name.s] = to_upper(seq->seq.s, seq->seq.l);
                ref_to_length[seq->name.s] = seq->seq.l;
            } 
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
                read_to_seq[seq->name.s] = to_upper(seq->seq.s, seq->seq.l);
                read_to_length[seq->name.s] = seq->seq.l;
            }
        }

        errtre << "Loaded " << read_to_seq.size() << " reads." << endl;

        cerr << errtre.str(); //<< "Loaded " << read_to_seq.size() << " reads." << endl;
        errtre.str("");
    }
    else{
        cerr << "Please provide a read file containing query sequences." << endl;
        exit(1);
    }

    if (kmer.size() < 1){
        #pragma omp master
        cerr << "No kmer size provided. The default of 16 will be used." << endl;
        kmer.push_back(16);
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

        #pragma omp for
        for (int i = 0; i < read_seq.size(); i++){
            tuple<hash_t*, int> hashes_and_num =  allhash_unsorted_64_fast(read_seq[i].second, read_to_length[read_seq[i].first], kmer);
            hash_t* hashes = std::get<0>(hashes_and_num);
            int num_hashes = std::get<1>(hashes_and_num);

            read_to_hashes[read_seq[i].first] = hashes;
            read_to_num_hashes[read_seq[i].first] = num_hashes;


            if (min_kmer_occ > 0){
                for (int j = 0; j < num_hashes; j++){
                    hash_t hashk = hashes[j];
                    #pragma omp atomic update
                    read_hash_to_depth[hashk]++;
                }
            }
        }

        if (max_samples < 1000000){
            // set ref_to_hashes to a new value
            // set ref_to_num_hashes to a new value
            // need a only_informative_kmers(hash_t*, int hash_num, int max_sample) function
        }

        #pragma omp for
        for (int i = 0; i < read_seq.size(); i++){
            // need a minhash_64_depth_filter(hash_t* hash, int hash_num, int sketch_size, map<read_hash_to_depth>, int maxDepth) function
            stringstream outre;

            hash_t* hashes = read_to_hashes[read_seq[i].first];
            int hash_len = read_to_num_hashes[read_seq[i].first];

            vector<hash_t> mins;
            if (min_kmer_occ > 0){

            }
            else{
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


    for (auto x : read_to_hashes){
        delete [] x.second;
    }
    for (auto y : ref_to_hashes){
        delete [] y.second;
    }

    // Make read and ref hashes
    //
    // make kmer to depth for reads
    //
    // make ref_only_informative kmers
    //
    // make minhash sketches
    //
    // compare and classify
    //
    // call



    kseq_destroy(seq);
    gzclose(fp);



    return 1;
}
