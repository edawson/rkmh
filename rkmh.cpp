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
    void print_help(char** argv){
        cerr << "Usage: " << argv[0] << " [options]" << endl
            << "Options:" << endl
            << "--reads/-r   <READFILE>" << endl
            << "--fasta/-f   <FASTAFILE>" << endl
            << "--kmer/-k    <KMERSIZE>" << endl
            << "--minhash/-m <HASHSIZE>" << endl
            << "--threads/-t <THREADS>" << endl
            << "--min-kmer-occurence/-D <MINOCCURENCE>" << endl
            << "--min-matches/-S <MINMATCHES>" << endl;
    }


int main(int argc, char** argv){
    char* ref_file;
    char* read_file;
    vector<char*> ref_files;
    vector<char*> read_files;
    vector<int> kmer;
    int sketch_size = -1;
    int threads = 1;
    int min_kmer_occ = 0;
    int min_matches = -1;

    stringstream errtre;
    map<string, string> ref_to_seq;
    map<string, vector<int64_t> > ref_to_hashes;
    map<string, vector<string> > ref_to_kmers;

    map<string, string> read_to_seq;
    map<string, vector<int64_t> > read_to_hashes;
    map<string, vector<string> > read_to_kmers;

    unordered_map<string, int> kmer_to_depth;
    unordered_map<int64_t, int> hash_to_depth;
    hash_to_depth.reserve(10000000);

    //int64_t d_arr [INT64_MAX]; array is too large


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
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "hm:k:r:f:t:D:S:", long_options, &option_index);
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

    if (ref_files.size() > 0){
        for (auto f : ref_files){
            fp = gzopen(f, "r");
            seq = kseq_init(fp);
            // Read in reads, cluster, spit it back out
            while ((l = kseq_read(seq)) >= 0) {
                ref_to_seq[seq->name.s] = to_upper(seq->seq.s);
                //ref_seq.push_back(std::make_pair(seq->name.s, to_upper(seq->seq.s)));
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
                read_to_seq[seq->name.s] = to_upper(seq->seq.s);
                //read_seq.push_back(make_pair(seq->name.s, seq->seq.s));
            }
        }

        errtre << "Loaded " << read_to_seq.size() << " reads." << endl;

        cerr << errtre.str(); //<< "Loaded " << read_to_seq.size() << " reads." << endl;
        errtre.str(std::string());
    }
    else{
        cerr << "Please provide a read file containing query sequences." << endl;
        exit(1);
    }

    if (kmer.size() < 1){
        cerr << "No kmer size provided. The default of 16 will be used." << endl;
        kmer.push_back(16);
    }

    vector<pair<string, string > > ref_seq(ref_to_seq.begin(), ref_to_seq.end());
    vector<pair<string, string> > read_seq(read_to_seq.begin(), read_to_seq.end());
    
    #pragma omp parallel
    {
    if (min_kmer_occ > 0){
        #pragma omp master
        cerr << "Making kmer depth map..." << endl;
        // fill in kmer_to_depth;
        //#pragma omp parallel for
        #pragma omp for
        for (int i = 0; i < read_seq.size(); i++){
            pair<string, string> n_to_s = read_seq[i];
            vector<int64_t> rhash = allhash_unsorted_64(n_to_s.second, kmer);
            //#pragma omp atomic
            #pragma omp critical
            {
                read_to_hashes[n_to_s.first] = rhash;
                for (int j = 0; j < rhash.size(); j++){
                    hash_to_depth[rhash[j]] += 1;
                }
            }
        }
    }


    if (sketch_size > 0){
        #pragma omp master
        cerr << "Making reference sketches..." << endl;

        //#pragma omp parallel for schedule(dynamic)
        #pragma omp for
        for (int i = 0; i < ref_seq.size(); i++){
            ref_to_hashes[ref_seq[i].first] = minhash_64(ref_seq[i].second, kmer, sketch_size, true);
        }

        #pragma omp master
        cerr << "Processed " << ref_to_hashes.size() << " references to MinHashes" << endl;

        //#pragma omp parallel for schedule(dynamic)
        #pragma omp for
        for (int i = 0; i < read_seq.size(); i++){
            stringstream outre;
            vector<int64_t> hashes;
            if (min_kmer_occ > 0){
                hashes = minhash_64_depth_filter(read_to_hashes[read_seq[i].first], sketch_size, true, min_kmer_occ, hash_to_depth);
            }
            else{
                hashes = minhash_64(read_seq[i].second, kmer, sketch_size, true);
            }
            tuple<string, int, int> result;
            result = classify_and_count(hashes, ref_to_hashes);
            
            
            bool depth_filter = hashes.size() == 0; 
            bool match_filter = std::get<1>(result) < min_matches;
            //#pragma omp critical
            outre  << "Sample: " << read_seq[i].first << "\t" << "Result: " << 
                std::get<0>(result) << "\t" << std::get<1>(result) << "\t" << std::get<2>(result) << "\t" <<
                (depth_filter ? "FAIL:DEPTH" : "") << "\t" << (match_filter ? "FAIL:MATCHES" : "") << endl;

            cout << outre.str();
        }
    }

    else{
        #pragma omp master
        cerr << "Performing direct kmer-based comparison." << endl;
        map<string, priority_queue<string>> readheap;
        map<string, string>::iterator itersk;
        //for (itersk = read_to_seq.begin(); itersk != read_to_seq.end(); itersk++){
            //ref_to_kmers[itersk->first] = multi_kmerize(itersk->second, kmer);
            //std::sort(ref_to_kmers[itersk->first].begin(), ref_to_kmers[itersk->first].end());
        //    readheap[itersk->first] = kmer_heap(itersk->second, kmer);
        //}
        map<string, priority_queue<string> > ref_heaps;
        for (itersk = ref_to_seq.begin(); itersk != ref_to_seq.end(); itersk++){
            ref_heaps[itersk->first] = kmer_heap(itersk->second, kmer);
        }
        
        errtre << "Processed " << ref_heaps.size() << " references to kmers." << endl;
        #pragma omp master
        cerr << errtre;

        //vector<pair<string, priority_queue<string> > > p_read_heaps(readheap.begin(), readheap.end());
        vector<pair<string, priority_queue<string> > > p_ref_heaps(ref_heaps.begin(), ref_heaps.end());
        #pragma omp for
        for (int i = 0; i < read_seq.size(); i++){

            stringstream outre;

            priority_queue<string> read_q = kmer_heap(read_seq[i].second, kmer);
            tuple<string, int, int> result = kmer_heap_classify(read_q, p_ref_heaps);
            outre << read_seq[i].first << "\t" << std::get<0>(result) << "\t" << std::get<1>(result) << "\t" << std::get<2>(result) << endl;
            cout << outre.str();
        }
        //errtre << "Processed " << ref_to_kmers.size() << " references to kmers.";
        //cerr << errtre.str() << endl; //<< "Processed " << ref_to_kmers.size() << " references to kmers." << endl;
/*
        for (itersk = read_to_seq.begin(); itersk != read_to_seq.end(); itersk++){
            read_to_kmers[itersk->first] = multi_kmerize(itersk->second, kmer);
            std::sort(read_to_kmers[itersk->first].begin(), read_to_kmers[itersk->first].end());

            tuple<string, int, int> result = kmer_classify(read_to_kmers[itersk->first], ref_to_kmers);
            cout  << "Sample: " << itersk->first << "\t"
                << "Result: " << std::get<0>(result) << "\t" << std::get<1>(result) << "\t" << std::get<2>(result) << endl;


        }
        //cerr << "Processed " << read_to_kmers.size() << " reads to kmers." << endl;
        */


    }
    }






    return 1;
}
