#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>
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
//#include "kseq.hpp"
#include "equiv.hpp"
#include "json.hpp"
#include "HASHTCounter.hpp"
#include "kseq_reader.hpp"

// for convenience
using json = nlohmann::json;

//KSEQ_INIT(gzFile, gzread)

    using namespace std;
    using namespace mkmh;
    using namespace HTC;
    using namespace KSR;



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
     * ./rkmh stream:
     *  sketch size: s
     *  kmer size: k
     *  threads: t
     *  fasta: f; behavior: hash these first and then wait for STDIN,
     *      hashing each fastq or string as it comes in (perhaps we should signal
     *      what's a fastq, what's fasta, and what's a sequence?)
     *  references: r
     *  min_kmer_occ, precomputed -M + -P
     *  minmatches: N
     *  min diff: D
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
     *
     * for calling, we might be interested in using some score
     *  to quantify how good a variant call is.
     *  We could use the number of adjacent supporting kmers
     *  and the depth / difference to average depth. We
     *  already do this to some degree by setting the depth
     *  threshold for a kmer to be considered a recovery at .9 * avg depth.
     *
     *  map<int, map<kmer, int> > pos_to_kmer_to_occurrence
     *  
     */

    vector<string> split(string s, char delim){
        vector<string> ret;
        stringstream sstream(s);
        string temp;
        while(getline(sstream, temp, delim)){
            ret.push_back(temp);
        }
        return ret;

    }

string join(vector<string> splits, string glue){                                              
    string ret = "";
    for (int i = 0; i < splits.size(); i++){
        if (i != 0){
            ret += glue;
        }
        ret += splits[i];
    }

    return ret;
}

void print_help(char** argv){
    cerr << "Usage: " << argv[0] << " { classify | call | hash | stream } [options]" << endl
        << "    classify: match each read to the reference it most closely resembles using MinHash sketches." << endl
        << "    call: determine the SNPs and 1-2bp INDELs that differ between a set of reads and their closest reference." << endl
        << "    hash: compute the MinHash sketches of a set of reads and/or references (for interop with Mash/sourmash)." << endl
        << "    stream: classify reads or sequences from STDIN. Low memory, real time, but possibly lower precision." << endl
        << "    filter: spit out reads that meet thresholds for match to ref, uniqueness, etc." << endl
        << endl;
}

void help_classify(char** argv){
    cerr << "Usage: " << argv[0] << " classify [options]" << endl
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

void help_call(char** argv){
    cerr << "Usage: " << argv[0] << " call [options]" << endl
        << "Options:" << endl
        << "--reference/-r <REF>      reference genomes in fasta format." << endl
        << "--fasta/-f <FASTA>        a fasta file to call mutations in relative to the reference." << endl
        << "--threads/-t <THREADS>    the number of OpenMP threads to utilize." << endl
        << "--window-len/-w <WINLEN>  the width of the sliding window to use for calculating average depth." << endl
        << "--depth/-d                output tab-separated values for position, avg depth, instantaneous depth, and rescued depth." << endl
        << endl;
}

void help_hash(char** argv){
    cerr << "Usage: " << argv[0] << " hash [options]" << endl
        << "Options:" << endl
        << "--fasta/-f  <FASTA>          fasta file to hash." << endl
        << "--reference/-r   <REF>       reference file to hash." << endl
        << "--sketch-size/-s <SKTCHSZ>   sketch size." << endl
        << "--kmer/-k <KMER>             kmer size to hash." << endl
        << "--min-kmer-occurrence <M>    Minimum kmer occurrence. Failing kmers are removed from sketch." << endl
        << "--min-informative/-I  <I>    Maximum number of samples a kmer can occur in before it is removed" << endl
        << "--threads/-t <THREADS>       number of OpenMP threads to utilize." << endl
        << "--wabbitize /-w              output Vowpal Wabbit compatible vectors" << endl;
}

void help_stream(char** argv){
    cerr << "Usage: " << argv[0] << " stream [options]" << endl
        << "Options:" << endl
        << "--input-stream/-i   classify reads coming from STDIN" << endl
        << "--output-reads/-z   output reads that pass the filters set, rather than classifications." << endl 
        << "--reference/-r   <REF>" << endl
        << "--fasta/-f   <FASTAFILE>" << endl
        << "--kmer/-k    <KMERSIZE>" << endl
        << "--sketch-size/-s <SKETCHSIZE>" << endl
        << "--ref-sketch / -S <REFSKTCHSZ>" << endl
        << "--threads/-t <THREADS>" << endl
        << "--min-kmer-occurence/-M <MINOCCURENCE>" << endl
        << "--min-matches/-N <MINMATCHES>" << endl
        << "--min-diff/-D    <MINDIFFERENCE>" << endl
        << "--min-informative/-I <MAXSAMPLES> only use kmers present in fewer than MAXSAMPLES" << endl
        << "--kmer-depth-map / -p <mapfile> the kmer depth map to use for min_kmer_occurence" << endl
        << "--ref-sample-map / -q <mapfile> the sample depth map for reference sample filtering." << endl
        << "--pre-fasta / -F  a file containing sketches in JSON format for reads." << endl
        << "--pre-reference / -R a file containing pre-hashed reference genomes in JSON format." << endl
        << endl;

}

void help_filter(char** argv){
    cerr << "Usage: " << argv[0] << " filter [options]" << endl
        << "Options: " << endl
        << "--input-stream/-i   classify reads coming from STDIN" << endl
        << "--output-reads/-z   output reads that pass the filters set, rather than classifications." << endl 
        << "--reference/-r   <REF>" << endl
        << "--fasta/-f   <FASTAFILE>" << endl
        << "--kmer/-k    <KMERSIZE>" << endl
        << "--sketch-size/-s <SKETCHSIZE>" << endl
        << "--ref-sketch / -S <REFSKTCHSZ>" << endl
        << "--threads/-t <THREADS>" << endl
        << "--min-kmer-occurence/-M <MINOCCURENCE>" << endl
        << "--min-matches/-N <MINMATCHES>" << endl
        << "--min-diff/-D    <MINDIFFERENCE>" << endl
        << "--min-informative/-I <MAXSAMPLES> only use kmers present in fewer than MAXSAMPLES" << endl
        << "--kmer-depth-map / -p <mapfile> the kmer depth map to use for min_kmer_occurence" << endl
        << "--ref-sample-map / -q <mapfile> the sample depth map for reference sample filtering." << endl
        << "--pre-fasta / -F  a file containing sketches in JSON format for reads." << endl
        << "--pre-reference / -R a file containing pre-hashed reference genomes in JSON format." << endl
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
            seq_keys.push_back( string(seq->name.s) );
            seq_seqs.push_back(x);
            seq_lens.push_back(seq->seq.l);
        } 
        gzclose(fp);
    }   
    kseq_destroy(seq);
}

void parse_fastas(vector<char*>& files,
        vector<string>& seq_keys,
        vector<char*>& seq_seqs,
        vector<int>& seq_lens,
        vector<string>& seq_quals){

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
            seq_quals.emplace_back(seq->qual.s);
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

#pragma omp parallel for
    for (int i = 0; i < keys.size(); i++){
        tuple<hash_t*, int> hashes_and_num =  allhash_unsorted_64_fast(name_to_seq[keys[i]], name_to_length[keys[i]], kmer);
        ret_to_hashes[keys[i]] = std::get<0>(hashes_and_num);
        ret_to_hash_num[keys[i]] = std::get<1>(hashes_and_num);
    }

}
void hash_sequences(vector<string>& keys,
        vector<char*>& seqs,
        vector<int>& lengths,
        vector<hash_t*>& hashes,
        vector<int>& hash_lengths,
        vector<int>& kmer,
        HASHTCounter& read_hash_counter,
        HASHTCounter& ref_hash_counter,
        bool doReadDepth,
        bool doReferenceDepth){


    if (doReadDepth){
#pragma omp parallel for
        for (int i = 0; i < keys.size(); i++){
            // Hash sequence
            tuple<hash_t*, int> hashes_and_num =  allhash_unsorted_64_fast(seqs[i], lengths[i], kmer);
            hashes[i] = std::get<0>(hashes_and_num);
            hash_lengths[i] = std::get<1>(hashes_and_num);
            // TODO this is awful. There has to be a safe way around it.
            //#pragma omp critical
            {
                for (int j = 0; j < hash_lengths[i]; j++){

                    //#pragma omp atomic update //TODO removing this is under testing
                    //++read_hash_counter[ (uint64_t) hashes[i][j] ];
                    read_hash_counter.increment( (uint64_t) hashes[i][j] );
                }
            }
        }
    }
    else if (doReferenceDepth){
#pragma omp parallel for
        for (int i = 0; i < keys.size(); i++){
            tuple<hash_t*, int> hashes_and_num =  allhash_unsorted_64_fast(seqs[i], lengths[i], kmer);
            hashes[i] = std::get<0>(hashes_and_num);
            hash_lengths[i] = std::get<1>(hashes_and_num);

            // create the set of hashes in the sample
            set<hash_t> sample_set (hashes[i], hashes[i] + hash_lengths[i]);
            //#pragma omp critical
            {
                for (auto x : sample_set){
                    //#pragma omp atomic update //TODO under testing
                    //++ref_hash_counter[ x ];
                    ref_hash_counter.increment( x );
                }
            }
        }
    }

    else{
#pragma omp parallel for
        for (int i = 0; i < keys.size(); i++){
            tuple<hash_t*, int> hashes_and_num =  allhash_unsorted_64_fast(seqs[i], lengths[i], kmer);
            hashes[i] = std::get<0>(hashes_and_num);
            hash_lengths[i] = std::get<1>(hashes_and_num);
        }

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

    if (doReadDepth){
#pragma omp parallel for
        for (int i = 0; i < keys.size(); i++){
            // Hash sequence
            tuple<hash_t*, int> hashes_and_num =  allhash_unsorted_64_fast(seqs[i], lengths[i], kmer);
            hashes[i] = std::get<0>(hashes_and_num);
            hash_lengths[i] = std::get<1>(hashes_and_num);
            // TODO this is awful. There has to be a safe way around it.
            {
                for (int j = 0; j < hash_lengths[i]; j++){
#pragma omp critical //TODO removing this is under testing
                    ++read_hash_to_depth[ hashes[i][j] ];
                }
            }
        }
    }
    else if (doReferenceDepth){
#pragma omp parallel for
        for (int i = 0; i < keys.size(); i++){
            tuple<hash_t*, int> hashes_and_num =  allhash_unsorted_64_fast(seqs[i], lengths[i], kmer);
            hashes[i] = std::get<0>(hashes_and_num);
            hash_lengths[i] = std::get<1>(hashes_and_num);

            // create the set of hashes in the sample
            set<hash_t> sample_set (hashes[i], hashes[i] + hash_lengths[i]);
            //#pragma omp critical
            {
                for (auto x : sample_set){
#pragma omp critical //TODO under testing
                    ++ref_to_sample_depth[x];
                }
            }
        }
    }

    else{
#pragma omp parallel for
        for (int i = 0; i < keys.size(); i++){
            tuple<hash_t*, int> hashes_and_num =  allhash_unsorted_64_fast(seqs[i], lengths[i], kmer);
            hashes[i] = std::get<0>(hashes_and_num);
            hash_lengths[i] = std::get<1>(hashes_and_num);
        }

    }

}

string sketch_to_json(string key,
                    vector<hash_t> mins,
                    int sketchlen,
                    vector<int> kmer,
                    int sketch_size){
    
}
void sketches_to_jsons(vector<string>& keys,
                        vector<vector<hash_t> >& mins,
                        vector<int>& sketchlens,
                        vector<int>& kmer,
                        int sketch_size){
    cout << "[" << endl;
    for (int i = 0; i < keys.size(); ++i){
        cout << sketch_to_json(keys[i], mins[i], sketchlens[i], kmer, sketch_size);
    }
    cout << "]" << endl;
}


void rkmh_binary_output(vector<string> keys,
                        vector<vector<hash_t> >& mins,
                        vector<int> sketchlens,
                        vector<int>& kmer,
                        int sketch_size){

}

void print_wabbit(string key,
                  vector<hash_t>& mins,
                  int sketch_len,
                  int sketch_size,
                  vector<int> counts,
                  vector<int> kmers,
                  string label = "XYX",
                  string nspace = "vir"){

       key = join( split(key, '|'), "_");
       cout << label << " 1.0 " << "`" << key << "|" << nspace;
       if (!counts.empty()){
            for (int i = 0; i < sketch_len; ++i){
                cout << " " << mins[i] << ":" << counts[i];
        }
       }
       else{
             for (int i = 0; i < sketch_len; ++i){
                cout << " " << mins[i] << ":1";
        }

        }
       cout << " |sketch k:" << kmers[0] << " " << "s:" << sketch_size;
              cout << endl;
}

json dump_hash_json(string key, int seqLen,
        vector<hash_t> mins,
        vector<int> kmer, 
        int sketchLen,
        string alphabet = "ATGC",
        string hash_type = "MurmurHash3_x64_128",
        bool canonical = true,
        int hash_bits = 64,
        int hash_seed = 42
        ){
    json j;
    j["name"] = key;

    stringstream kstr;
    for (int i = 0; i < kmer.size(); i++){
        kstr << kmer[i];
        if (i < kmer.size() - 1){
            kstr << " ";
        }
    }
    j["kmer"] = kstr.str();
    j["alphabet"] = alphabet;
    j["preserveCase"] = "false";
    j["canonical"] = (canonical ? "true" : "false");
    j["hashType"] = hash_type;
    j["hashBits"] = hash_bits;
    j["hashSeed"] = hash_seed;
    j["seqLen"] = seqLen;
    j["sketches"] = {
        {"name", key},
        {"length", sketchLen},
        {"comment", ""},
        {"hashes", mins}
    };

    return j;
}

json dump_hashes(vector<string> keys,
        vector<int> seqlens,
        vector<vector<hash_t>> hashes,
        vector<int> kmer,
        int sketch_size){

    json j;
    for (int i = 0; i < keys.size(); i++){
        j.push_back({
            {"name", keys[i]},
            {"alphabet", "ATGC"},
            {"canonical", "false"},
            {"hashBits", 64},
            {"hash_type", "MurmurHash3_x64_128"},
            {"hash_seed", 42},
            {"sketches", hashes[i]},
            {"length", sketch_size},
            {"kmer", kmer},
            {"preserveCase", "false"}
        });
    }

    return j;
}

std::tuple<vector<string> ,
    vector<int> ,
    vector<vector<hash_t> > ,
    vector<int> ,
    int> load_hashes(json jj){
        cout << jj << endl;
        cerr << "Loading not implemented" << endl;
        exit(1);

    }

std::tuple<vector<string>,
    vector<int> ,
    vector<vector<hash_t> >,
    vector<int> ,
    int> load_hashes(istream ifi){
        json jj;
        ifi >> jj;
        return load_hashes(jj);
    }




std::tuple<vector<string> ,
    vector<int> ,
    vector<vector<hash_t> >,
    vector<int>,
    int> load_hashes(string filename){

    }

int main_stream(int argc, char** argv){
    vector<char*> ref_files;
    vector<char*> read_files;
    vector<char*> pre_read_files;
    vector<char*> pre_ref_files;

    vector<int> kmer;

    int sketch_size = 1000;
    int threads = 1;
    int min_kmer_occ = -1;
    int min_matches = -1;
    int min_diff = 0;
    int max_samples = 100000;

    string read_kmer_map_file = "";
    string ref_kmer_map_file = "";

    bool doReadDepth = false;
    bool doReferenceDepth = false;

    bool useHASHTs = false;
    int ref_sketch_size = 0;

    bool streamify_me_capn = false;
    bool output_reads = false;
    bool merge_sketch = false;

    // TODO still need:
    // prehashed depth map for reads/ref
    // prehashed reads / refs
    //

    int c;
    int optind = 2;

    if (argc <= 2){
        help_stream(argv);
        exit(1);
    }

    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"kmer", no_argument, 0, 'k'},
            {"fasta", required_argument, 0, 'f'},
            {"reference", required_argument, 0, 'r'},
            {"sketch-size", required_argument, 0, 's'},
            {"ref-sketch", required_argument, 0, 'S'},
            {"threads", required_argument, 0, 't'},
            {"min-kmer-occurence", required_argument, 0, 'M'},
            {"min-matches", required_argument, 0, 'N'},
            {"min-diff", required_argument, 0, 'D'},
            {"max-samples", required_argument, 0, 'I'},
            {"pre-reads", required_argument, 0, 'F'},
            {"pre-references", required_argument, 0, 'R'},
            {"read-kmer-map-file", required_argument, 0, 'p'},
            {"ref-kmer-map-file", required_argument, 0, 'q'},
            {"in-stream", no_argument, 0, 'i'},
            {"output-reads", no_argument, 0, 'z'},
            {"merge-sketch", no_argument, 0, 'm'},
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "zmhdk:f:r:s:S:t:M:N:I:R:F:p:q:iD:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case 'm':
                merge_sketch = true;
                break;
            case 'F':
                pre_read_files.push_back(optarg);
                break;
            case 'R':
                pre_ref_files.push_back(optarg);
                break;
            case 'p':
                read_kmer_map_file = optarg;
                break;
            case 'q':
                ref_kmer_map_file = optarg;
                break;
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
            case 'N':
                min_matches = atoi(optarg);
                break;
            case 'D':
                min_diff = atoi(optarg);
                break;
            case '?':
            case 'h':
                print_help(argv);
                exit(1);
                break;
            case 's':
                sketch_size = atoi(optarg);
                break;
            case 'S':
                useHASHTs = true;
                ref_sketch_size = 3 * atoi(optarg);
                break;
            case 'M':
                min_kmer_occ = atoi(optarg);
                doReadDepth = true;
                break;
            case 'I':
                max_samples = atoi(optarg);
                doReferenceDepth = true;
                break;
            case 'i':
                streamify_me_capn = true;
                break;
            case 'z':
                output_reads = true;
                break;
            default:
                print_help(argv);
                abort();

        }
    }

    if (sketch_size == -1){
        cerr << "Sketch size unset." << endl
            << "Will use the default sketch size of s = 1000" << endl;
        sketch_size = 1000;
    }

    if (kmer.size() == 0){
        cerr << "No kmer size(s) provided. Will use a default kmer size of 16." << endl;
        kmer.push_back(16);
    }


    omp_set_num_threads(threads);
    // Read in depth map for reads and refs if provided
    HASHTCounter read_hash_counter(200000000);
    HASHTCounter ref_hash_counter(200000000);
    if (!read_kmer_map_file.empty()){
        ifstream ifi(read_kmer_map_file);
        string xline;
        int xsz = 0;
        int pos = 0;
        while (getline(ifi, xline, '\n')){
            if (pos != 0){
                read_hash_counter.set(pos, stoi(xline));
            }
            else{
                xsz = stoi(xline);
                read_hash_counter.size(xsz);
            }
        }
    }

    if (!ref_kmer_map_file.empty()){

    }
    //read in prehashed sequences
    if (!pre_read_files.empty()){

    }
    if (!pre_ref_files.empty()){

    }

    vector<string> ref_keys;
    vector<char*> ref_seqs;
    vector<int> ref_lens;

    vector<string> read_keys;
    vector<char*> read_seqs;
    vector<int> read_lens;


    //read in new refs/reads, hash them and keep their sketches.

    if (!ref_files.empty()){
        parse_fastas(ref_files, ref_keys, ref_seqs, ref_lens);
    }
    if (!read_files.empty()){
        parse_fastas(read_files, read_keys, read_seqs, read_lens);
    }
    else if (read_files.empty() && (!streamify_me_capn)){
        cerr << "No reads provided. Please provide a read file with -f or pipe reads in with -i" << endl;
        exit(1);
    }

    vector<hash_t*> ref_hashes(ref_keys.size());
    vector<int> ref_hash_lens(ref_keys.size());

    vector<hash_t*> read_hashes(read_keys.size());
    vector<int> read_hash_lens(read_keys.size());

    vector<hash_t*> ref_mins(ref_keys.size());
    int* ref_min_lens = new int [ref_keys.size()];
    int* ref_min_starts = new int [ref_keys.size()];

    vector<hash_t*> read_mins(read_keys.size());
    int* read_min_starts = new int [ read_keys.size() ];
    int* read_min_lens = new int [read_keys.size() ];

    vector<vector<hash_t> > mergeable(read_keys.size());

    int ref_htc_size = 0;
    int read_htc_size = 0;

    // Attempt to make the HTC load 40% or less
    // Otherwise we'll get collisions and that's no good
    if (!read_keys.empty() && read_kmer_map_file.empty()){
        read_htc_size = (int)  ((double) (read_keys.size() * read_lens[0]) * 5);

        //read_hash_counter.size( read_htc_size );
    }
    if (!ref_keys.empty() && ref_kmer_map_file.empty()){
        ref_htc_size = (int) ((double) (ref_keys.size() * ref_lens[0]) * 5);

        //ref_hash_counter.size( ref_htc_size );
    }

    vector<vector<string> > results(threads);

    if (!ref_files.empty()){
        hash_sequences(ref_keys, ref_seqs, ref_lens, ref_hashes, ref_hash_lens, kmer, read_hash_counter, ref_hash_counter, false, doReferenceDepth);
    }


    if (!read_files.empty()){
        hash_sequences(read_keys, read_seqs, read_lens, read_hashes, read_hash_lens, kmer, read_hash_counter, ref_hash_counter, doReadDepth, false);
    }

    //Time to calculate mins for references!
    int nthreads = 0;
    int tid = 0;
#pragma omp parallel private(tid, nthreads)
    {
        tid = omp_get_thread_num();
#pragma omp for
        for (int i = 0; i < ref_keys.size(); i++){
            ref_min_starts[i] = 0;
            ref_min_lens[i] = 0;
            ref_mins[i] = new hash_t[ sketch_size ];
            std::sort(ref_hashes[i], ref_hashes[i] + ref_hash_lens[i]);
            if (max_samples < 100000){
                for (int j = 0; j < ref_hash_lens[i], ref_min_lens[i] < sketch_size; ++j){
                    //if (ref_sketch_lens[i] >= sketch_size){
                    //    break;
                    //}
                    hash_t curr = *(ref_hashes[i] + j);
                    //cerr << ref_hash_counter.get(curr) << endl;
                    if (curr != 0 && ref_hash_counter.get(curr) <= max_samples){
                        ref_mins[i][ref_min_lens[i]] = curr;
                        //ref_mins[i][ref_min_lens[i]] = ref_hashes[i][j];
                        ++ref_min_lens[i];
                        if (ref_min_lens[i] == sketch_size){
                            break;
                        }
                    }
                    else{
                        continue;
                    }
                }

            }
            else{
                while (ref_hashes[i][ref_min_starts[i]] == 0 && ref_min_starts[i] < ref_hash_lens[i]){
                    ++ref_min_starts[i];
                }
                for (int j = ref_min_starts[i]; j < ref_hash_lens[i], ref_min_lens[i] < sketch_size; ++j){
                    *(ref_mins[i] +ref_min_lens[i]) = *(ref_hashes[i] + j);
                    ++ref_min_lens[i];
                }

            }
            ref_min_starts[i] = 0;

            delete [] ref_hashes[i];

        }

        // Classify existing reads
        // conveniently, read_keys.size() will be zero if there are no reads.
#pragma omp for
        for (int i = 0; i < read_keys.size(); i++){
            stringstream outre;
            read_mins[i] = new hash_t[ sketch_size ];
            read_min_lens[i] = 0;
            read_min_starts[i] = 0;
            std::sort(read_hashes[i], read_hashes[i] + read_hash_lens[i]);
            if (doReadDepth){
                for (int j = 0; j < read_hash_lens[i]; ++j){

                    if (read_hashes[i][j] != 0 && read_hash_counter.get(read_hashes[i][j]) >= min_kmer_occ){
                        read_mins[i][read_min_lens[i]] = *(read_hashes[i] + j);

                        ++(read_min_lens[i]);
                        if (read_min_lens[i] == sketch_size){
                            break;
                        }
                    }
                    else{
                        continue;
                    }
                }
            }
            else{
                while (read_hashes[i][ read_min_starts[i] ] == 0 && read_min_starts[i] < read_hash_lens[i]){
                    ++read_min_starts[i];
                }
                for (int j = read_min_starts[i]; j < read_hash_lens[i]; ++j){
                    read_mins[i][read_min_lens[i]] = *(read_hashes[i] + j);
                    ++(read_min_lens[i]);
                    if (read_min_lens[i] == sketch_size){
                        break;
                    }
                }
            }
            read_min_starts[i] = 0;
            delete [] read_hashes[i];

            if (!merge_sketch){
                tuple<string, int, int, bool> result;
                result = classify_and_count_diff_filter(ref_keys, ref_mins, read_mins[i], ref_min_starts, read_min_starts[i], ref_min_lens, read_min_lens[i], sketch_size, min_diff);


                bool depth_filter = read_min_lens[i] <= 0; 
                bool match_filter = std::get<1>(result) < min_matches;

                outre  << "Sample: " << read_keys[i] << "\t" << "Result: " << 
                    std::get<0>(result) << "\t" << std::get<1>(result) << "\t" << std::get<2>(result) << "\t" <<
                    (depth_filter ? "FAIL:DEPTH" : "") << "\t" << (match_filter ? "FAIL:MATCHES" : "") << "\t" << (std::get<3>(result) ? "" : "FAIL:DIFF") << endl;

            //#pragma omp critical
            //cout << outre.str();
                tid = omp_get_thread_num();
                results[tid].push_back(outre.str());
                outre.str("");
                delete [] read_mins[i];
                }
                else{
                    vector<hash_t> read_vec(read_mins[i], read_mins[i] + read_min_lens[i]);
                    mergeable[i] = read_vec;
                }


        }

        if (!merge_sketch){
            #pragma omp single
            {
                for (int i = 0; i < results.size(); ++i){
                    for (int j = 0; j < results[i].size(); ++j){
                        cout << results[i][j];
                    }
                }
            }

        }
        else{
#pragma omp single
            {
            tuple<hash_t*, int> merged = merge(mergeable, sketch_size);
            tuple<string, int, int, bool> result;
            result = classify_and_count_diff_filter(ref_keys, ref_mins, std::get<0>(merged), ref_min_starts, 0, ref_min_lens, std::get<1>(merged), sketch_size, min_diff);


            bool depth_filter = std::get<1>(merged) <= 0; 
            bool match_filter = std::get<1>(result) < min_matches;
            stringstream outre;
            outre  << "Sample: " << "merged" << "\t" << "Result: " << 
            std::get<0>(result) << "\t" << std::get<1>(result) << "\t" << std::get<2>(result) << "\t" <<
            (depth_filter ? "FAIL:DEPTH" : "") << "\t" << (match_filter ? "FAIL:MATCHES" : "") << "\t" << (std::get<3>(result) ? "" : "FAIL:DIFF") << endl;
            cout << outre.str();
            }
        }

        // Take in a quartet of lines from STDIN (FASTQ format??)
        // or perhaps just individual read sequences and names (or give them names dynamically
        // hash them, possibly with kmer depth filtering,
        // create a sketch, then classify and report the classification.
        //
        // Thanks heavens for https://biowize.wordpress.com/2013/03/05/using-kseq-h-with-stdin/

        if (streamify_me_capn){
            FILE *instream = NULL;
            instream = stdin;

            gzFile fp = gzdopen(fileno(instream), "r");
            kseq_t *seq = kseq_init(fp);
            stringstream outre;
            while (kseq_read(seq) >= 0){
                to_upper(seq->seq.s, seq->seq.l);

                string name = string(seq->name.s);
                int len = seq->seq.l;

                // hash me
                tuple<hash_t*, int> hashes_and_num =  allhash_unsorted_64_fast(seq->seq.s, len, kmer);
                stringstream outre;

                hash_t* hashes = std::get<0>(hashes_and_num);
                int hashlen = std::get<1>(hashes_and_num);

                std::sort(hashes, hashes + hashlen);
                // and then just sketch me
                // TODO need to handle some read_depth
                int sketch_start = 0;
                int sketch_len = 0;
                hash_t* mins = new hash_t[sketch_size];
                if (min_kmer_occ > 0){
                    for (int i = 0; i < hashlen; ++i){
                        hash_t curr = *(hashes + i);
                        if (read_hash_counter.get(curr) >= min_kmer_occ && curr != 0){
                            mins[sketch_len] = curr;
                            ++sketch_len;
                        }
                        if (sketch_len == sketch_size){
                            break;
                        }
                    }
                }
                else{
                    while (hashes[sketch_start] == 0 && sketch_start < hashlen){
                        ++sketch_start;
                    }
                    for (int i = sketch_start; i < hashlen; ++i){
                        mins[sketch_len++] = *(hashes + i);
                        if (sketch_len == sketch_size){
                            break;
                        }
                    }
                }
                sketch_start = 0;
                // so I can get my
                // classification
                tuple<string, int, int, bool> result;
                result = classify_and_count_diff_filter(ref_keys, ref_mins, mins, ref_min_starts, sketch_start, ref_min_lens, sketch_len, sketch_size, min_diff);

                bool depth_filter = sketch_len <= 0; 
                bool match_filter = std::get<1>(result) < min_matches;

                if (output_reads){
                    if (!depth_filter && !match_filter && std::get<3>(result)){
                        outre << ">" << seq->name.s << endl
                            << seq->seq.s << endl
                            << "+" << endl
                            << seq->qual.s << endl;
                    }
                }
                else{
                    outre  << "Sample: " << name << "\t" << "Result: " << 
                        std::get<0>(result) << "\t" << std::get<1>(result) << "\t" << std::get<2>(result) << "\t" <<
                        (depth_filter ? "FAIL:DEPTH" : "") << "\t" << (match_filter ? "FAIL:MATCHES" : "") << "\t" << (std::get<3>(result) ? "" : "FAIL:DIFF") << endl;
                }
#pragma omp critical
                cout << outre.str();
                outre.str("");

                delete [] hashes;
                delete [] mins;

                // classification
            }
            kseq_destroy(seq);
            gzclose(fp);

        }
    }


    delete [] ref_min_lens;
    delete [] read_min_lens;
    delete [] read_min_starts;
    delete [] ref_min_starts;

    return 0;



}


/**
 * Filter indidvidual sequencing reads based on their % match to a panel
 * of reference MinHash sketches.
 * */
int main_filter(int argc, char** argv){
    vector<char*> ref_files;
    vector<char*> read_files;
    vector<char*> pre_read_files;
    vector<char*> pre_ref_files;

    vector<int> kmer;

    int sketch_size = 1000;
    int threads = 1;
    int min_kmer_occ = -1;
    int min_matches = -1;
    int min_diff = 0;
    int max_samples = 100000;

    string read_kmer_map_file = "";
    string ref_kmer_map_file = "";

    bool doReadDepth = false;
    bool doReferenceDepth = false;

    bool useHASHTs = false;
    int ref_sketch_size = 0;

    bool streamify_me_capn = false;
    bool output_reads = false;

    // TODO still need:
    // prehashed depth map for reads/ref
    // prehashed reads / refs
    //

    int c;
    int optind = 2;

    if (argc <= 2){
        help_filter(argv);
        exit(1);
    }

    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"kmer", no_argument, 0, 'k'},
            {"fasta", required_argument, 0, 'f'},
            {"reference", required_argument, 0, 'r'},
            {"sketch-size", required_argument, 0, 's'},
            {"ref-sketch", required_argument, 0, 'S'},
            {"threads", required_argument, 0, 't'},
            {"min-kmer-occurence", required_argument, 0, 'M'},
            {"min-matches", required_argument, 0, 'N'},
            {"min-diff", required_argument, 0, 'D'},
            {"max-samples", required_argument, 0, 'I'},
            {"pre-reads", required_argument, 0, 'F'},
            {"pre-references", required_argument, 0, 'R'},
            {"read-kmer-map-file", required_argument, 0, 'p'},
            {"ref-kmer-map-file", required_argument, 0, 'q'},
            {"in-stream", no_argument, 0, 'i'},
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "hdk:f:r:s:S:t:M:N:I:R:F:p:q:iD:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case 'F':
                pre_read_files.push_back(optarg);
                break;
            case 'R':
                pre_ref_files.push_back(optarg);
                break;
            case 'p':
                read_kmer_map_file = optarg;
                break;
            case 'q':
                ref_kmer_map_file = optarg;
                break;
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
            case 'N':
                min_matches = atoi(optarg);
                break;
            case 'D':
                min_diff = atoi(optarg);
                break;
            case '?':
            case 'h':
                print_help(argv);
                exit(1);
                break;
            case 's':
                sketch_size = atoi(optarg);
                break;
            case 'S':
                useHASHTs = true;
                ref_sketch_size = 3 * atoi(optarg);
                break;
            case 'M':
                min_kmer_occ = atoi(optarg);
                doReadDepth = true;
                break;
            case 'I':
                max_samples = atoi(optarg);
                doReferenceDepth = true;
                break;
            case 'i':
                streamify_me_capn = true;
                break;
            default:
                print_help(argv);
                abort();

        }
    }

    if (sketch_size == -1){
        cerr << "Sketch size unset." << endl
            << "Will use the default sketch size of s = 10000" << endl;
        sketch_size = 1000;
    }

    if (kmer.size() == 0){
        cerr << "No kmer size(s) provided. Will use a default kmer size of 16." << endl;
        kmer.push_back(16);
    }


    omp_set_num_threads(threads);
    // Read in depth map for reads and refs if provided
    if (!read_kmer_map_file.empty()){

    }
    if (!ref_kmer_map_file.empty()){

    }
    //read in prehashed sequences
    if (!pre_read_files.empty()){

    }
    if (!pre_ref_files.empty()){

    }

    vector<string> ref_keys;
    vector<char*> ref_seqs;
    vector<int> ref_lens;

    vector<string> read_keys;
    vector<char*> read_seqs;
    vector<string> read_quals;
    vector<int> read_lens;


    //read in new refs/reads, hash them and keep their sketches.

    if (!ref_files.empty()){
        parse_fastas(ref_files, ref_keys, ref_seqs, ref_lens);
    }
    if (!read_files.empty()){
        parse_fastas(read_files, read_keys, read_seqs, read_lens, read_quals);
    }

    vector<hash_t*> ref_hashes(ref_keys.size());
    vector<int> ref_hash_lens(ref_keys.size());

    vector<hash_t*> read_hashes(read_keys.size());
    vector<int> read_hash_lens(read_keys.size());

    vector<hash_t*> ref_mins(ref_keys.size());
    int* ref_min_lens = new int [ref_keys.size()];
    int* ref_min_starts = new int [ref_keys.size()];

    vector<hash_t*> read_mins(read_keys.size());
    int* read_min_starts = new int [ read_keys.size() ];
    int* read_min_lens = new int [read_keys.size() ];


    HASHTCounter read_hash_counter(10000000);
    HASHTCounter ref_hash_counter(10000000);

    vector<vector<string> > results(threads);

    if (!ref_files.empty()){
        hash_sequences(ref_keys, ref_seqs, ref_lens, ref_hashes, ref_hash_lens, kmer, read_hash_counter, ref_hash_counter, false, doReferenceDepth);
    }


    if (!read_files.empty()){
        hash_sequences(read_keys, read_seqs, read_lens, read_hashes, read_hash_lens, kmer, read_hash_counter, ref_hash_counter, doReadDepth, false);
    }

    //Time to calculate mins for references!
#pragma omp parallel
    {
        //int tid = omp_get_thread_num();
#pragma omp for
        for (int i = 0; i < ref_keys.size(); i++){
            ref_min_starts[i] = 0;
            ref_min_lens[i] = 0;
            ref_mins[i] = new hash_t[ sketch_size ];
            std::sort(ref_hashes[i], ref_hashes[i] + ref_hash_lens[i]);
            if (max_samples < 100000){
                for (int j = 0; j < ref_hash_lens[i], ref_min_lens[i] < sketch_size; ++j){
                    //if (ref_sketch_lens[i] >= sketch_size){
                    //    break;
                    //}
                    hash_t curr = *(ref_hashes[i] + j);
                    //cerr << ref_hash_counter.get(curr) << endl;
                    if (curr != 0 && ref_hash_counter.get(curr) <= max_samples){
                        ref_mins[i][ref_min_lens[i]] = curr;
                        //ref_mins[i][ref_min_lens[i]] = ref_hashes[i][j];
                        ++ref_min_lens[i];
                        if (ref_min_lens[i] == sketch_size){
                            break;
                        }
                    }
                    else{
                        continue;
                    }
                }

            }
            else{
                while (ref_hashes[i][ref_min_starts[i]] == 0 && ref_min_starts[i] < ref_hash_lens[i]){
                    ++ref_min_starts[i];
                }
                for (int j = ref_min_starts[i]; j < ref_hash_lens[i], ref_min_lens[i] < sketch_size; ++j){
                    *(ref_mins[i] +ref_min_lens[i]) = *(ref_hashes[i] + j);
                    ++ref_min_lens[i];
                }

            }
            ref_min_starts[i] = 0;

            delete [] ref_hashes[i];

        }

        // Classify existing reads
        // conveniently, read_keys.size() will be zero if there are no reads.
#pragma omp for
        for (int i = 0; i < read_keys.size(); i++){
            stringstream outre;
            read_mins[i] = new hash_t[ sketch_size ];
            read_min_lens[i] = 0;
            read_min_starts[i] = 0;
            std::sort(read_hashes[i], read_hashes[i] + read_hash_lens[i]);
            if (doReadDepth){
                for (int j = 0; j < read_hash_lens[i]; ++j){

                    if (read_hashes[i][j] != 0 && read_hash_counter.get(read_hashes[i][j]) >= min_kmer_occ){
                        read_mins[i][read_min_lens[i]] = *(read_hashes[i] + j);

                        ++(read_min_lens[i]);
                        if (read_min_lens[i] == sketch_size){
                            break;
                        }
                    }
                    else{
                        continue;
                    }
                }
            }
            else{
                while (read_hashes[i][ read_min_starts[i] ] == 0 && read_min_starts[i] < read_hash_lens[i]){
                    ++read_min_starts[i];
                }
                for (int j = read_min_starts[i]; j < read_hash_lens[i]; ++j){
                    read_mins[i][read_min_lens[i]] = *(read_hashes[i] + j);
                    ++(read_min_lens[i]);
                    if (read_min_lens[i] == sketch_size){
                        break;
                    }
                }
            }
            read_min_starts[i] = 0;
            delete [] read_hashes[i];

            tuple<string, int, int, bool> result;
            result = classify_and_count_diff_filter(ref_keys, ref_mins, read_mins[i], ref_min_starts, read_min_starts[i], ref_min_lens, read_min_lens[i], sketch_size, min_diff);


            bool depth_filter = read_min_lens[i] <= 0; 
            bool match_filter = std::get<1>(result) < min_matches;

            //cerr << read_keys[i] << " " << read_seqs[i] << endl
            //    << read_quals[i] << endl;

            if (!depth_filter && !match_filter && std::get<3>(result)){
                outre << ">" << string(read_keys[i]) << endl
                    << string(read_seqs[i], read_lens[i]) << endl
                    << "+" << endl
                    << string(read_quals[i]) << endl;

                //results[tid].push_back(outre.str());
                #pragma omp critical
                {
                    cout << outre.str();
                    outre.str("");
                }
            }
            delete [] read_mins[i];


        }
    }
    /**for (int i = 0; i < results.size(); ++i){
      for (int j = 0; j < results[i].size(); ++j){
      cout << results[i][j];
      }
      }**/

    // Take in a quartet of lines from STDIN (FASTQ format??)
    // or perhaps just individual read sequences and names (or give them names dynamically
    // hash them, possibly with kmer depth filtering,
    // create a sketch, then classify and report the classification.
    //
    // Thanks heavens for https://biowize.wordpress.com/2013/03/05/using-kseq-h-with-stdin/

    if (streamify_me_capn){
        FILE *instream = NULL;
        instream = stdin;

        gzFile fp = gzdopen(fileno(instream), "r");
        kseq_t *seq = kseq_init(fp);
        stringstream outre;
        while (kseq_read(seq) >= 0){
            to_upper(seq->seq.s, seq->seq.l);

            string name = string(seq->name.s);
            int len = seq->seq.l;

            // hash me
            tuple<hash_t*, int> hashes_and_num =  allhash_unsorted_64_fast(seq->seq.s, len, kmer);
            stringstream outre;

            hash_t* hashes = std::get<0>(hashes_and_num);
            int hashlen = std::get<1>(hashes_and_num);

            std::sort(hashes, hashes + hashlen);
            // and then just sketch me
            // TODO need to handle some read_depth
            int sketch_start = 0;
            int sketch_len = 0;
            hash_t* mins = new hash_t[sketch_size];
            if (min_kmer_occ > 0){
                for (int i = 0; i < hashlen; ++i){
                    hash_t curr = *(hashes + i);
                    if (read_hash_counter.get(curr) >= min_kmer_occ && curr != 0){
                        mins[sketch_len] = curr;
                        ++sketch_len;
                    }
                    if (sketch_len == sketch_size){
                        break;
                    }
                }
            }
            else{
                while (hashes[sketch_start] == 0 && sketch_start < hashlen){
                    ++sketch_start;
                }
                for (int i = sketch_start; i < hashlen; ++i){
                    mins[sketch_len++] = *(hashes + i);
                    if (sketch_len == sketch_size){
                        break;
                    }
                }
            }
            sketch_start = 0;
            // so I can get my
            // classification
            tuple<string, int, int, bool> result;
            result = classify_and_count_diff_filter(ref_keys, ref_mins, mins, ref_min_starts, sketch_start, ref_min_lens, sketch_len, sketch_size, min_diff);

            bool depth_filter = sketch_len <= 0; 
            bool match_filter = std::get<1>(result) < min_matches;

            if (!depth_filter && !match_filter && std::get<3>(result)){
                outre << ">" << seq->name.s << endl
                    << seq->seq.s << endl
                    << "+" << endl
                    << seq->qual.s << endl;
            }
#pragma omp critical
            {
                cout << outre.str();
                outre.str("");
            }
            delete [] hashes;
            delete [] mins;

            // classification
        }
        kseq_destroy(seq);
        gzclose(fp);
        if (streamify_me_capn){
            //#pragma omp parallel
            //{
#pragma omp single nowait
            {
                FILE *instream = NULL;
                instream = stdin;

                gzFile fp = gzdopen(fileno(instream), "r");
                kseq_t *seq = kseq_init(fp);


                while (kseq_read(seq) >= 0){

                    to_upper(seq->seq.s, seq->seq.l);

                    string name = string(seq->name.s);
                    int len = seq->seq.l;

                    // hash me

#pragma omp task_async
                    {
                        tuple<hash_t*, int> hashes_and_num =  allhash_unsorted_64_fast(seq->seq.s, len, kmer);
                        stringstream outre;

                        hash_t* hashes = std::get<0>(hashes_and_num);
                        int hashlen = std::get<1>(hashes_and_num);

                        std::sort(hashes, hashes + hashlen);
                        // and then just sketch me
                        // TODO need to handle some read_depth
                        int sketch_start = 0;
                        int sketch_len = 0;
                        hash_t* mins = new hash_t[sketch_size];
                        if (min_kmer_occ > 0){
                            for (int i = 0; i < hashlen; ++i){
                                hash_t curr = *(hashes + i);
                                if (read_hash_counter.get(curr) >= min_kmer_occ && curr != 0){
                                    mins[sketch_len] = curr;
                                    ++sketch_len;
                                }
                                if (sketch_len == sketch_size){
                                    break;
                                }
                            }
                        }
                        else{
                            while (hashes[sketch_start] == 0 && sketch_start < hashlen){
                                ++sketch_start;
                            }
                            for (int i = sketch_start; i < hashlen; ++i){
                                mins[sketch_len++] = *(hashes + i);
                                if (sketch_len == sketch_size){
                                    break;
                                }
                            }
                        }
                        sketch_start = 0;
                        // so I can get my
                        // classification
                        tuple<string, int, int, bool> result;
                        result = classify_and_count_diff_filter(ref_keys, ref_mins, mins, ref_min_starts, sketch_start, ref_min_lens, sketch_len, sketch_size, min_diff);

                        bool depth_filter = sketch_len <= 0; 
                        bool match_filter = std::get<1>(result) < min_matches;

                        outre  << "Sample: " << name << "\t" << "Result: " << 
                            std::get<0>(result) << "\t" << std::get<1>(result) << "\t" << std::get<2>(result) << "\t" <<
                            (depth_filter ? "FAIL:DEPTH" : "") << "\t" << (match_filter ? "FAIL:MATCHES" : "") << "\t" << (std::get<3>(result) ? "" : "FAIL:DIFF") << endl;

                        //#pragma omp critical
                        cout << outre.str();
                        outre.str("");

                        delete [] hashes;
                        delete [] mins;
                    }
                }

                kseq_destroy(seq);
                gzclose(fp);
            }
        }
        }


        delete [] ref_min_lens;
        delete [] read_min_lens;
        delete [] read_min_starts;
        delete [] ref_min_starts;



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

        vector<int> kmer;

        int sketch_size = 1000;
        int threads = 1;
        int min_kmer_occ = 1;
        int min_matches = -1;
        int min_diff = 0;
        int max_samples = 1000000;
        int window_len = 100;

        bool show_depth = false;
        bool output_vcf = true;

        int c;
        int optind = 2;

        if (argc <= 2){
            help_call(argv);
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
                {"show-depth", required_argument, 0, 'd'},
                {"max-samples", required_argument, 0, 'I'},
                {"window-len", required_argument, 0, 'w'},
                {0,0,0,0}
            };

            int option_index = 0;
            c = getopt_long(argc, argv, "hdk:f:r:s:t:M:N:I:w:", long_options, &option_index);
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
                case 'I':
                    max_samples = atoi(optarg);
                    break;
                case 'w':
                    window_len = atoi(optarg);
                    break;
                case 'd':
                    show_depth = true;
                    output_vcf = false;
                    break;
                default:
                    print_help(argv);
                    abort();

            }
        }

        if (sketch_size == -1){
            cerr << "Sketch size unset." << endl
                << "Will use the default sketch size of s = 1000" << endl;
            sketch_size = 1000;
        }

        if (kmer.size() == 0){
            cerr << "No kmer size(s) provided. Will use a default kmer size of 16." << endl;
            kmer.push_back(16);
        }
        else if (kmer.size() > 1){
            cerr << "Only a single kmer size may be used for calling." << endl
                << "Sizes provided: ";
            for (auto k : kmer){
                cerr << k << " ";
            }
            cerr << endl;
            cerr << "Please choose a single kmer size." << endl;
            exit(1);
        }

        omp_set_num_threads(threads);

        //TODO switch to c arrays from vectors?
        // we know the size of these and we carry the lengths around

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

        // TODO try something faster than a map / unordered_map?
        // a massive array?

        unordered_map<hash_t, int> read_hash_to_depth;
        read_hash_to_depth.reserve(1000000);
        unordered_map<hash_t, int> ref_hash_to_num_samples;
        ref_hash_to_num_samples.reserve(1000000);


#pragma omp master
        cerr << "Parsing sequences..." << endl;

        if (ref_files.size() >= 1){
            parse_fastas(ref_files, ref_keys, ref_seqs, ref_lens);
        }
        else{
            cerr << "No references were provided. Please provide at least one reference file in fasta/fastq format." << endl;
            help_classify(argv);
            exit(1);
        }

        if (read_files.size() >= 1){
            parse_fastas(read_files, read_keys, read_seqs, read_lens);
        }
        else{
            cerr << "No reads were provided. Please provide at least one read file in fasta/fastq format." << endl;
            help_classify(argv);
            exit(1);
        }



#pragma omp master
        cerr << " Done." << endl <<
            ref_keys.size() << " references and " << read_keys.size() << " reads parsed." << endl;

        vector<hash_t*> ref_hashes(ref_keys.size());
        vector<int> ref_hash_nums(ref_keys.size());

        vector<hash_t*> read_hashes(read_keys.size());
        vector<int> read_hash_nums(read_keys.size());


        vector<vector<hash_t> > ref_mins(ref_keys.size(), vector<hash_t>(1));


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

        std::function<double(vector<int>)> avg = [](vector<int> n_list){
            int ret = 0;
            for (int x = 0; x < n_list.size(); x++){
                ret += n_list.at(x);
            }
            return (double) ret / (double) n_list.size();
        };

        std::function<double(int* depth_arr, int start, int window, int size)> avg_arr = [](int* depth_arr, int start, int window, int size){
            int end =  1;
            return 2.0;
        };

        vector<char> a_ret = {'C', 'T', 'G'};
        vector<char> c_ret = {'T', 'G', 'A'};
        vector<char> t_ret = {'C', 'G', 'A'};
        vector<char> g_ret = {'A', 'C', 'T'};


        std::function<vector<char>(char)> rotate_snps = [&a_ret, &c_ret, &g_ret, &t_ret](char c){

            if ( c == 'A' || c == 'a'){
                return a_ret;
            }
            else if (c == 'T' || c == 't'){
                return t_ret;
            }
            else if (c == 'C' || c == 'c'){
                return c_ret;
            }
            else if (c == 'G' || c == 'g'){
                return g_ret; 
            }

        };

        std::function<vector<string>(string)> permute = [&](string x){
            //int k_pos;
            //string mut;

            vector<string> ret;
            // SNPs, 3 * sequence length
            for (int i = 0; i < x.size(); i++){
                char orig = x[i];
                vector<char> other_chars = rotate_snps(x[i]);
                for (int j = 0; j < other_chars.size(); j++){
                    x[i] = other_chars[j];
                    ret.push_back(x);
                }
                x[i] = orig;
            }

            //DELs, 1bp, 1 * sequence length
            // TODO probably doesn't delete the first base correctly
            for (int i = 0; i < x.size() - 1; i++){
                char orig = x[i];
                stringstream tmp;
                for (int strpos = 0; strpos < x.size(); strpos++){
                    if (strpos != i){
                        tmp << x[strpos];
                    }
                }
                ret.push_back(tmp.str());
            }

            //DELs, 2bp, 1 * sequence length
            // TODO probably doesn't work on first occurence correctly either
            for (int i = 0; i < x.size() - 2; i++){
                char orig = x[i];
                stringstream tmp;
                for (int strpos = 0; strpos < x.size(); strpos++){
                    if (strpos != i & strpos != i + 1){
                        tmp << x[strpos];
                    }
                }
                ret.push_back(tmp.str());
            }

            // 1bp insertions, 1 * sequence length
            for (int i = 0; i < x.size(); i++){
                char orig = x[i];
                stringstream tmp;
                for (int strpos = 0; strpos < x.size(); strpos++){
                    tmp << x[strpos]; 
                }
            }

            // 2bp insertions, 1 * sequence length
            for (int i = 0; i < x.size(); i++){
                char orig = x[i];
                stringstream tmp;
                for (int strpos = 0; strpos < x.size(); strpos++){

                }
            }

            // Homopolymer extension / contraction?

            return ret;
        };

        /***
         * TODO: try using the entire sequence as the window size,
         * calculating average depth. Maybe store depth windows rather than calculate them?
         *
         * Also reverse complement entire sequence and use the (seqlen - i) trick to compare the same kmers,
         * only calling reverse_reverse_complement once per sequence
         *
         * map[ reference ] -> vector<depth>: windowed depth at each position
         * vector<array<int> > 
         *  conveniently the same length as the reference.
         */


        vector<int*> depth_arrs(ref_keys.size());

        if (ref_keys.size() > 1){
            cerr << "WARNING: more than one ref provided. VCF will not be correct" << endl;
            
        }
        if (output_vcf){
            cout << "##fileformat=VCF4.2\n##source=rkmh\n##reference=" << ref_files[0] << endl <<
                "##INFO=<ID=KD,Number=1,Type=Integer,Description=\"Number of times call for specific kmer appears\">" << endl
                << "##INFO=<ID=MD,Number=1,Type=Integer,Description=\"Maximum depth found for the rescue kmer.\">" << endl
                << "##INFO=<ID=RD,Number=1,Type=Integer,Description=\"Average depth in region\">"
                << "##INFO=<ID=OD,Number=1,Type=Integer,Description=\"Depth of original kmer at site before modification.\">"
                << endl;
        }

        std::function<string(vector<string>)> join = [](vector<string> strs){
            stringstream strstream;
            for (int i = 0; i < strs.size() - 1; i++){
                strstream << strs[i] << "_";
            }
            strstream << strs[strs.size() - 1];
            return strstream.str();
        };

        map<string, int> call_count;
        map<string, int> call_max_depth;
        map<string, int> call_avg_depth;
        map<string, int> call_orig_depth;
        vector<string> outbuf;
        outbuf.reserve(1000);


#pragma omp parallel
        {

            list<int> d_window;


            #pragma omp for
            for (int i = 0; i < ref_keys.size(); i++){
                // This loop iterates over the reference genomes.

                stringstream outre;
                //outre << ref_keys[i] << endl;
                //cout << outre.str(); outre.str("");

                for (int j = 0; j < ref_hash_nums[i]; j++){
                    // This loop iterates over a single reference genome
                    // i.e. its sequence

                    // This is a hacky way of calculating depth using a sliding window.
                    int depth = read_hash_to_depth[ref_hashes[i][j]];
                    d_window.push_back(depth);
                    if (d_window.size() > window_len){
                        d_window.pop_front();
                    }

                    int avg_d = avg(vector<int>(d_window.begin(), d_window.end()));
                    int max_rescue = 0;

                    //int max_rescue = 0;

                    // This line outputs the current avg depth at a position.
                    if (show_depth){
                        outre << j << "\t" << avg_d << "\t" <<  depth;

                    }
                    if (depth < .5 * avg_d){
                        string ref = string(ref_seqs[i] + j, kmer[0]);
                        string alt(ref);
                        string d_alt = (j > 0) ? string(ref_seqs[i] + j - 1, kmer[0] + 1) : "";

                        // SNPs
                        for (int alt_pos = 0; alt_pos < alt.size(); alt_pos++){
                            char orig = alt[alt_pos];
                            for (auto x : rotate_snps(orig)){
                                alt[alt_pos] = x;
                                int alt_depth = read_hash_to_depth[calc_hash(alt)];
                                max_rescue = max_rescue > alt_depth ? max_rescue : alt_depth;

                                if ( !show_depth && alt_depth >= .1 * avg_d & alt_depth > depth){
                                    int pos = j + alt_pos + 1;
                                    //outre << "CALL: " << orig << "->" << x << "\t" << "POS: " << pos << "\tRESCUE_DEPTH: " << alt_depth << "\tORIGINAL_DEPTH: " << depth << "\tSURROUNDING_AVG: " << avg_d<< endl;
                                    if (output_vcf){
                                        #pragma omp critical
                                        {
                                            stringstream sstream;
                                            sstream << ref_keys[i] << "\t" << pos << "\t" <<
                                                "." << "\t" << orig << "\t" << x;
                                            string s = sstream.str();
                                            call_count[s] += 1;
                                            call_avg_depth[s] = max(avg_d, call_avg_depth[s]);
                                            call_orig_depth[s] = max(call_orig_depth[s], depth);
                                            if (alt_depth > call_max_depth[s]){
                                                call_max_depth[s] = alt_depth;
                                            }
                                        }
                                    }
                                    else{
                                        outre << "CALL: " << orig << "->" << x << "\t" << "POS: " << pos << "\tRESCUE_DEPTH: " << alt_depth << endl;
                                        outre << "\t" << "old: " << ref << endl << "\t" << "new: " << alt;
                                    }
                                }
                                alt[alt_pos] = orig;
                            }

                        }

                        char atgc[4] = {'A', 'T', 'G', 'C'};

                        // Deletions
                        // TODO both insertions and deletions are tough because we don't know whether to take a trailing or
                        // precending character in building the new kmer. Either might be optimal.
                        if (j > 0){
                            for (int alt_pos = 1; alt_pos < d_alt.size(); alt_pos++){
                                stringstream mod;
                                char orig = d_alt[alt_pos];
                                mod << d_alt.substr(0, alt_pos) << d_alt.substr(alt_pos + 1, d_alt.length() - alt_pos);
                                int alt_depth = read_hash_to_depth[calc_hash(mod.str())];
                                if (output_vcf && alt_depth > 0.9 * avg_d){
                                    int pos = j + alt_pos + 1;
                                    stringstream sstream;
                                    sstream << ref_keys[i] << "\t" << pos << "\t" << "." << "\t" << orig << "\t" << "-";
                                    string s = sstream.str();
                                    call_count[s] += 1;
                                    call_avg_depth[s] = max(call_avg_depth[s], avg_d);
                                    call_orig_depth[s] = max(call_orig_depth[s], depth);
                                    if (alt_depth > call_max_depth[s]){
                                        call_max_depth[s] = alt_depth;
                                    }
                                }
                            }
    
                            // Insertions
                            //
                        }
                        


                    }

                    if (show_depth){
                        outre << "\t" << (max_rescue > 0 ? max_rescue : depth);
                    }

                }

            }

        }

        for (auto x : call_count){
            cout << x.first << "\t" << "99" << "\t" << "PASS" << "\t" << "KC=" << x.second << ";" <<
                "MD=" << call_max_depth[x.first] << ";" << "RD=" << call_avg_depth[x.first] << ";OD=" << call_orig_depth[x.first] << endl;
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
     * TODO: streaming kmer counting
     */
    int main_hash(int argc, char** argv){
        vector<char*> ref_files;
        vector<char*> read_files;


        vector<int> kmer;

        int sketch_size = 0;
        int threads = 1;
        int min_kmer_occ = 0;
        int min_matches = -1;
        int min_diff = 0;
        int max_samples = 100000;

        bool doReadDepth = false;
        bool doReferenceDepth = false;
        bool wabbitize = false;
        bool merge_sketch = false;

        bool traditional_minhash = false;
        bool output_kmers = false;
        bool output_counts = false;
        string outname = "";

        int c;
        int optind = 2;

        if (argc <= 2){
            help_hash(argv);
            exit(1);
        }

        while (true){
            static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"count", no_argument, 0, 'c'},
                {"kmer", no_argument, 0, 'k'},
                {"output=kmers", no_argument, 0, 'K'},
                {"merge-sample", no_argument, 0, 'm'},
                {"wabbitize", no_argument, 0, 'w'},
                {"fasta", required_argument, 0, 'f'},
                {"reference", required_argument, 0, 'r'},
                {"sketch-size", required_argument, 0, 's'},
                {"threads", required_argument, 0, 't'},
                {"min-kmer-occurence", required_argument, 0, 'M'},
                {"max-samples", required_argument, 0, 'I'},
                {"out-prefix", required_argument, 0, 'o'},
                {0,0,0,0}
            };

            int option_index = 0;

            c = getopt_long(argc, argv, "ThcwKk:f:r:s:t:mM:I:o:", long_options, &option_index);
            if (c == -1){
                break;
            }

            switch (c){
                case 'T':
                    traditional_minhash = true;
                    break;
                case 'c':
                    output_counts = true;
                    break;
                case 'm':
                    merge_sketch = true;
                    break;
                case 'w':
                    wabbitize = true;
                    break;
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
                case 'K':
                    output_kmers = true;
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
                    doReadDepth = true;
                    break;
                case 'I':
                    max_samples = atoi(optarg);
                    doReferenceDepth = true;
                    break;
                case 'o':
                    outname = string(optarg);
                    break;
                default:
                    print_help(argv);
                    abort();

            }
        }

        if (kmer.size() != 0){
            cerr << "Using a kmer size of " << kmer[0] << endl; 
        }
        else {
            cerr << "Using default kmer size of 16." << endl;
            kmer.push_back(16);
        }

        omp_set_num_threads(threads);

        // In parallel
        // read in a chunk of fastq reads
        // hash those reads
        // output said reads

        
        HASHTCounter read_counter(6400000000);
        kseq_t *seq;
        char* f = read_files[0];
        gzFile fp;
        int l;
        fp = gzopen(f, "r");
        seq = kseq_init(fp);
        while ((l = kseq_read(seq)) >= 0) {
            to_upper(seq->seq.s, seq->seq.l);
            char * x = new char[seq->seq.l];
            memcpy(x, seq->seq.s, seq->seq.l);
            int len = seq->seq.l;
            cout << seq->name.s << "\t";
            if (!output_kmers){
                tuple<hash_t*, int> r = allhash_unsorted_64_fast(x, len, kmer);
                
                int nhash = std::get<1>(r);
                hash_t* hashes = std::get<0>(r);
                for (int j = 0; j < nhash; ++j){
                    cout << "\t" << *(hashes + j);
                }
                cout << endl;
            }
            else{
                // mkmh_kmer_list_t kl = kmerize(x, len, kmer[0]);
                // for (int i = 0; i < kl.length; ++i){
                //     cout << "\t" << *(kl.kmers + i);
                // }
                // cout << endl;
                print_kmers(x, len, kmer[0]);
                cout << endl;
            }
            
            delete [] x;
        } 
        gzclose(fp);
        kseq_destroy(seq);
        return 0;
    }


    /**
     *  Search: take in a reference set of hashes and query
     *  a set of FASTA/FASTQ files against them to find queries
     *  that contain those hashes.
     * 
     */
    int main_search(int argc, char** argv){
        // These are just lists of hashes,
        // one list per line
        // token[0] is the name of the reference
        vector<char*> ref_files;
        // These are fastqs/fastas
        vector<char*> read_files;

        int bufsz = 500;

        string summary_file_suffix = "rkmh.summary.txt";

        int threads = 1;
        vector<int> kmer;

        int c;
        int optind = 2;

        if (argc <= 2){
            help_hash(argv);
            exit(1);
        }

        while (true){
            static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"fasta", required_argument, 0, 'f'},
                {"reference", required_argument, 0, 'r'},
                {"threads", required_argument, 0, 't'},
                {"kmer", required_argument, 0, 'k'},
                {0,0,0,0}
            };

            int option_index = 0;

            c = getopt_long(argc, argv, "k:f:r:h", long_options, &option_index);
            if (c == -1){
                break;
            }

            switch (c){
                case 't':
                    threads = atoi(optarg);
                    break;
                case 'k':
                    kmer.push_back(stoi(optarg));
                    break;
                case 'f':
                    read_files.push_back(optarg);
                    break;
                case 'r':
                    ref_files.push_back(optarg);
                    break;
                case '?':
                case 'h':
                    //print_help(argv);
                    exit(1);
                    break;
                default:
                    //print_help(argv);
                    abort();

            }
        }

        std::map<string, HASHTCounter*> ref_to_ht;
        for (auto r : ref_files){
            ifstream infile(r);
            string line;
            while(getline(infile, line)){
                vector<string> tokens = split(line, '\t');
                HASHTCounter* refhtc = new HASHTCounter(65000);
                ref_to_ht[tokens[0]] = refhtc;
                for (int i = 2; i < tokens.size(); ++i){
                    ref_to_ht[tokens[0]]->increment((hash_t) stoull(tokens[i]));
                }
            }
        }

        

        for (auto f : read_files){
            KSEQ_Reader kh;
            kh.open(f);
            kh.buffer_size(bufsz);
            ksequence_t* buf;
            int buflen;
            int l = kh.get_next_buff(buf, buflen);
            while (l == 0){
                // Iterate over a buffer of sequences
                for (int i = 0; i < buflen; ++i){
                    stringstream seqstr;
                    seqstr << (buf + i)->name << "\t";
                    int seqlen = (buf + i)->length;
                    tuple<hash_t*, int> hashes = allhash_unsorted_64_fast((buf + i)->sequence, seqlen, kmer);
                    hash_t* start = std::get<0>(hashes);
                    int num = std::get<1>(hashes);
                    for (map<string, HASHTCounter*>::iterator it = ref_to_ht.begin();
                    it != ref_to_ht.end();
                    it++){
                        int intersection = 0;
                        for (int j = 0; j < num; j++){
                            if (it->second->get( *(start + j) ) > 0){
                                intersection++;
                            }
                        }
                        cout << it->first << ":" << intersection << "/" << "NA" << endl;
                    }
                }
             l= kh.get_next_buff(buf, buflen);           
            }
        }
        return 0;
    }



    /**
     * Counts kmers and outputs their counts in a map
     * that can be used for streaming.
     * Essentially, it's streaming main_hash
     * main_count( int argc, char** argv){
     *  
     * }
     */
    int main_count(int argc, char** argv){
        vector<char*> read_files;
        vector<int> kmer;
        int threads = 1;

        int bz = 1000;

        int c;
        int optind = 2;

        if (argc <= 2){
            help_hash(argv);
            exit(1);
        }

        while (true){
            static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"fasta", required_argument, 0, 'f'},
                {"threads", required_argument, 0, 't'},
                {0,0,0,0}
            };

            int option_index = 0;

            c = getopt_long(argc, argv, "k:tf:h", long_options, &option_index);
            if (c == -1){
                break;
            }

            switch (c){
                case 't':
                    threads = atoi(optarg);
                    break;
                case 'k':
                    kmer.push_back(stoi(optarg));
                    break;
                case 'f':
                    read_files.push_back(optarg);
                    break;
                case '?':
                case 'h':
                    //print_help(argv);
                    exit(1);
                    break;
                default:
                    //print_help(argv);
                    abort();

            }
        }

        omp_set_num_threads(threads);
        HASHTCounter htc(640000);
        //KSEQ_Reader* kt = new KSEQ_Reader();
        //kt->buffer_size(bz);

    /**#pragma omp parallel
    {
    #pragma omp single
    {
                for (auto r : read_files){
                    //int l = 0;
                    //kt->open(r);
                    //while (l == 0){
                    //{
                        //vector<string> seqs;
                        //l = kt->get_next_sequence_buffer(seqs);
                            //for (int i = 0; i < seqs.size(); ++i){
                                    //#pragma omp task shared(seqs)
                                    //{
                                    //string s = seqs[i];
                                    //int nums = seqs[i].size();
                                    //tuple<hash_t*, int> r = allhash_unsorted_64_fast(s.c_str(), nums, kmer);
                                    //hash_t* h = std::get<0>(r);
                                    //int n = std::get<1>(r);
                                    //for (int i = 0; i < n; ++i){
                                        //htc.increment( *(h + i) );
                                    //}
                                    //}
                                //}
                            //}

                    //}
                }
            }
    }
    **/
    //delete kt;

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
            help_classify(argv);
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

        if (sketch_size == -1){
            cerr << "Sketch size unset." << endl
                << "Will use the default sketch size of s = 1000" << endl;
            sketch_size = 1000;
        }

        if (kmer.size() == 0){
            cerr << "No kmer size(s) provided. Will use a default kmer size of 16." << endl;
            kmer.push_back(16);
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

        unordered_map<hash_t, int> read_hash_to_depth;
        read_hash_to_depth.reserve(10000);
        unordered_map<hash_t, int> ref_hash_to_num_samples;
        ref_hash_to_num_samples.reserve(10000);


#pragma omp master
        cerr << "Parsing sequences..." << endl;

        if (ref_files.size() >= 1){
            parse_fastas(ref_files, ref_keys, ref_seqs, ref_lens);
        }
        else{
            cerr << "No references were provided. Please provide at least one reference file in fasta/fastq format." << endl;
            help_classify(argv);
            exit(1);
        }

        if (read_files.size() >= 1){
            parse_fastas(read_files, read_keys, read_seqs, read_lens);
        }
        else{
            cerr << "No reads were provided. Please provide at least one read file in fasta/fastq format." << endl;
            help_classify(argv);
            exit(1);
        }



#pragma omp master
        cerr << " Done." << endl <<
            ref_keys.size() << " references and " << read_keys.size() << " reads parsed." << endl;

        vector<hash_t*> ref_hashes(ref_keys.size());
        vector<int> ref_hash_nums(ref_keys.size());

        vector<hash_t*> read_hashes(read_keys.size());
        vector<int> read_hash_nums(read_keys.size());


        int* read_sketch_starts = new int[(read_keys.size())];
        int* ref_sketch_starts = new int [(ref_keys.size())];
        int* read_sketch_lens = new int [(read_keys.size())];
        int* ref_sketch_lens = new int [(ref_keys.size())];
        vector<hash_t*> read_sketches(read_keys.size());
        vector<hash_t*> ref_sketches(ref_keys.size());

        int sbuf_size = 100;

        vector<string> sbuf(sbuf_size);
        int sbuf_ind = 0;

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

#pragma omp parallel
        {

#pragma omp for

            for (int i = 0; i < ref_keys.size(); i++){
                ref_sketch_starts[i] = 0;
                ref_sketch_lens[i] = 0;
                ref_sketches[i] = new hash_t[ sketch_size ];
                std::sort(ref_hashes[i], ref_hashes[i] + ref_hash_nums[i]);
                if (max_samples < 10000){
                    for (int j = 0; j < ref_hash_nums[i]; j++){
                        if (ref_hashes[i][j] != 0 && ref_hash_to_num_samples[ref_hashes[i][j]] <= max_samples)
                            ref_sketches[i][ref_sketch_lens[i]] = ref_hashes[i][j];
                        ++ref_sketch_lens[i];

                        if (ref_sketch_lens[i] >= sketch_size){
                            break;
                        }
                        else{
                            continue;
                        }
                    }

                }
                else{
                    ref_sketches[i] = ref_hashes[i];
                    while (ref_hashes[i][ref_sketch_starts[i]] == 0 && ref_sketch_starts[i] < ref_hash_nums[i]){
                        ++ref_sketch_starts[i];
                    }
                    int diff = ref_hash_nums[i] - ref_sketch_starts[i];
                    ref_sketch_lens[i] = (diff < sketch_size) ? diff : sketch_size; 

                }

            }

#pragma omp for
            for (int i = 0; i < read_keys.size(); i++){
                read_sketch_starts[i] = 0;
                read_sketch_lens[i] = 0;
                hash_t* hh = read_hashes[i];
                int hh_l = read_hash_nums[i];
                stringstream outre;

                sort(hh, hh + hh_l);
                if (min_kmer_occ > 0){
                    read_sketches[i] = new hash_t[sketch_size];
                    while(hh[read_sketch_starts[i]] == 0 && read_sketch_starts[i] < read_hash_nums[i]){
                        ++(read_sketch_starts[i]);
                    }

                    for (int d_ind = read_sketch_starts[i]; d_ind < read_hash_nums[i]; ++d_ind){
                        int k_depth = 0;
                        // I'm sorry for using auto here, but
                        // read_hash_to_depth is a map<hash_t, int> FWIW.
                        auto k_d_iter = read_hash_to_depth.find(hh[d_ind]);
                        if (k_d_iter != read_hash_to_depth.end()){
                            k_depth = k_d_iter->second;
                        }
                        if ( k_depth >= min_kmer_occ){
                            read_sketches[i][read_sketch_lens[i]] = hh[d_ind];
                            ++(read_sketch_lens[i]);
                            if (read_sketch_lens[i] == sketch_size){
                                break;
                            }
                        }
                    }
                    read_sketch_starts[i] = 0;

                }
                else{
                    read_sketches[i] = read_hashes[i];
                    while(read_sketches[i][read_sketch_starts[i]] == 0 && read_sketch_starts[i] < read_hash_nums[i]){
                        ++read_sketch_starts[i];
                    }
                    int diff = hh_l - read_sketch_starts[i];
                    read_sketch_lens[i] = diff < sketch_size ? diff : sketch_size;

                }
                tuple<string, int, int, bool> result;
                result = classify_and_count_diff_filter(ref_keys, ref_sketches, read_sketches[i], ref_sketch_starts, read_sketch_starts[i], ref_sketch_lens, read_sketch_lens[i], sketch_size, min_diff);


                bool depth_filter = read_sketch_lens[i] <= 0; 
                bool match_filter = std::get<1>(result) < min_matches;

                outre  << "Sample: " << read_keys[i] << "\t" << "Result: " << 
                    std::get<0>(result) << "\t" << std::get<1>(result) << "\t" << std::get<2>(result) << "\t" <<
                    (depth_filter ? "FAIL:DEPTH" : "") << "\t" << (match_filter ? "FAIL:MATCHES" : "") << "\t" << (std::get<3>(result) ? "" : "FAIL:DIFF") << endl;




                #pragma omp critical
                cout << outre.str();

                outre.str("");
                if (read_hashes[i] != read_sketches[i]){
                    delete [] read_hashes[i];
                }
                delete [] read_sketches[i];


            }

            }

            for (int sbi = 0; sbi < sbuf_ind; ++sbi){
                cout << sbuf[sbi];
            }


            for (auto x : read_seqs){
                delete [] x;
            }
            for (int i = 0; i < ref_keys.size(); i++){
                if (ref_hashes[i] != ref_sketches[i]){
                    delete [] ref_hashes[i];
                }
                delete [] ref_sketches[i];
            }

            for (auto y : ref_seqs){
                delete [] y;
            }

            delete [] read_sketch_lens;
            delete [] ref_sketch_lens;

            delete [] read_sketch_starts;
            delete [] ref_sketch_starts;
            return 0;
        }


        int main(int argc, char** argv){

            if (argc <= 1){
                print_help(argv);
                exit(1);
            }
            string cmd = argv[1];
            if (cmd == "classify"){
                return main_classify(argc, argv);
            }
            else if (cmd == "hash"){
                return main_hash(argc, argv);
            }
            else if (cmd == "call"){
                return main_call(argc, argv);
            }
            else if (cmd == "count"){
                return main_count(argc, argv);
            }
            else if (cmd == "search"){
                return main_search(argc, argv);
            }
            else if (cmd == "stream"){
                return main_stream(argc, argv);
            }
            else if (cmd == "filter"){
                return main_filter(argc, argv);
            }
            else{
                print_help(argv);
                exit(1);
            }

        }
