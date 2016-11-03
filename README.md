RKMH: Read classification by Kmers / MinHash
--------------------------------------------
Eric T Dawson  
June 2016


[![](https://images.microbadger.com/badges/image/erictdawson/rkmh.svg)](https://microbadger.com/images/erictdawson/rkmh "Get your own image badge on microbadger.com")

### What is it
rkmh performs identification of *individual long reads* and alignment-free variant calling
using MinHash (as implemented in [Mash](https://github.com/marbl/Mash)).


We're using rkmh to identify which strains are present in infections with multiple strains of the same virus and any mutations they might have.
rkmh could also be used to mark reads as possible contaminants or call mutations in novel strains that have been
sequenced. You could also cluster the output of `rkmh hash` to discover separations between sequences in your sample.

### License
MIT, but please cite the repository if you use it.

### Dependencies and build process
The only external dependencies should be zlib and a compiler supporting OpenMP. To download and build:  

                    git clone --recursive https://github.com/edawson/rkmh.git  
                    cd rkmh  
                    make  

This should build rkmh and its library dependencies (mkmh and murmur3).


### *NEW*: Stream
rkmh can now stream reads through, using roughly constant memory. This is useful in read filtration workflows.
This command functions almost identically to `classify` and performs the same read classification task.

The `-M` flag will have no effect (currently) on stream reads. I'm working on methods to prehash the read depth filters
and pass them as files to rkmh.

`-I` still works, however references are still built in memory. We advise not doing big reference (e.g. the Human genome)
quite yet. Wait for the next release.

### Classify 
rkmh requires a set of reads and a set of references in the FASTA/FASTQ format. Reads need not
be in FASTQ format.


To use MinHash sketch of size 1000, and a kmer size of 10:  
```./rkmh classify -r references.fa -f reads.fq -k 10 -s 1000```

There's also now a filter for minimum kmer occurrence in a read set, compatible with the MinHash sketch.
To only use kmers that occur more than 10 times in the reads:  
```./rkmh classify -r references.fa -f reads.fq -k 10 -s 1000 -M 100```

There is also a filter that will fail reads with fewer than some number of matches to any reference.
It's availble via the `-N` flag:  
```./rkmh -r references.fa -f reads.fq -k 10 -s 1000 -M 100 -N 10```


**A note on optimum kmer size**: we've had a lot of success with k <= 15 on data fron ONT's R7 pore. I don't have any R9 flowcells around lab, but 
I expect we'll do a bit better on R9 given what others have been showing off.

### Call
Once you've identified which reference a set of reads most closely matches, you may want to figure out the differences between your set of reads
and your reference. `rkmh call` uses a brute-force approach to produce a list of candidate mutations / sequencing errors present in a readset.

```rkmh call -r ref.fa -f reads.fq -k 12 -t 4```  

We advise using only one reference during call, as it's relatively slow (~10x longer than classification, 10 seconds for 1100 reads). For example, you might first classify your reads using `classify`, then
for the top classification in your set run `rkmh call`.

### Hash
You might want to see the hashes generated by rkmh for debugging purposes. To do so, use the `hash` command.

```rkmh hash -r ref.fa -f reads.fq -k 12 -s 1000``` 


### Other options
These are extra options for the `classify` and `hash` commands. Some of them are also applicable to `call`. For full usage, just
type `./rkmh` or `./rkmh <command>` at the command line to get the help message.


```-t / --threads <INT>               number of OpenMP threads to use (default is 1)```  
```-M / --min-kmer-occurence <INT>    minimum number of times a kmer must appear in the set of reads to be included in a read's MinHash sketch.```  
```-N / --min-matches <INT>           minimum number of matches a read must have to any reference to be considered classified.```  
```-I / --max-samples <INT>           remove kmers that appear in more than <INT> reference genomes.```  
```-D / --min-difference <INT>        flag reads that have two matches within <INT> hashes of each other as failing.```   
```-k / --kmer <INT>                  the kmer size to use for hashing. Multiple kmer sizes may be passed, but they must all use the -k <INT> format (i.e. -k 12 -k 14 -k 16...)```   
```-s / --sketch-size                 the number of hashes to use when comparing reads / references.```    
```-f / --fasta                       a FASTA/FASTQ file to use as a read set. Can be passed multiple times (i.e. -f first.fa -f second.fa...)``` 
```-r / --reference                   a FASTA/FASTQ file to use as a reference set. Can be passed multiple times (i.e. -r ref.fa -r ref_second.fa...)```   



### Performance
On a set of 1000 minION reads from a known HPV strain, rkmh is ~97% accurate (correctly placing the read in the right strain
of 182 input reference strains) and runs in <20 seconds. With the kmer depth and minimum match filters we're approaching 100% accuracy for about the same run time.
We're working on ways to improve sensitivity with further filtering and correction.


rkmh is threaded using OpenMP. Hashing can handle more than 400 reads/second (400 * 7kb means we're running over 2,500,000 basepairs / second), with some room still left for improvement.


We've tested up to 100,000 6.5kb reads + 182 references in a bit over 8GB of RAM, but we're working to scale to larger genomes and more reads. We've run an E. coli
run (actually, Nick Loman's R7.3 dataset against 6 E. coli references) on a desktop with 16GB of RAM. We think with a few tweaks we can do a lot better.


### Getting help
Please post to the [github](https://github.com/edawson/rkmh.git) for help.
