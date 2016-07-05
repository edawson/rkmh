RKMH: Read classification by Kmers or MinHash
--------------------------------------------
Eric T Dawson  
June 2016

### What is it
rkmh performs MinHash (as implemented in [Mash](https://github.com/marbl/Mash)) on a set of reads
and reference genomes to classify each read in the set to its nearest match in the reference set.
It can also classify reads by directly comparing all of their kmers.

### License
MIT, but please cite the repository if you use it.

### Dependencies and build process
The only external dependencies should be zlib and a compiler supporting OpenMP. To download and build:  

                    git clone --recursive https://github.com/edawson/rkmh.git  
                    cd rkmh  
                    make  

This should build rkmh and its library dependencies (mkmh and murmur3).

### Usage
rkmh required a set of reads and a set of references in the FASTA/FASTQ format. Reads need not
be in FASTQ format.

To classify reads by kmers of length 10:  
```./rkmh -f references.fa -r reads.fq -k 10```

To do the same, but use a MinHash sketch of size 1000 instead of just comparing kmers:  
```./rkmh -f references.fa -r reads.fq -k 10 -m 1000```

There's also now a filter for minimum kmer occurrence in a read set, compatible with the MinHash sketch.
To only use kmers that occur more than 10 times in the reads:  
```./rkmh -f references.fa -r reads.fq -k 10 -m 1000 -D 100```

There is also a filter that will fail reads with fewer than some number of matches to any reference.
It's availble via the `-P` flag:  
```./rkmh -f references.fa -r reads.fq -k 10 -m 1000 -D 100 -S 10```

### Other options
```-t / --threads               number of OpenMP threads to use (default is 1)```  
```-D / --min-kmer-occurence    minimum number of times a kmer must appear in the set of reads to be included in a read's MinHash sketch.```  
```-S / --min-matches           minimum number of matches a read must have to any reference to be considered classified.```  

### Performance
On a set of 1000 minION reads from a known HPV strain, rkmh is ~97% accurate (correctly placing the read in the right strain
of 182 input reference strains) and runs in <30 seconds. With the kmer depth and minimum match filters we're approaching >99% accuracy for about the same run time.
We're working on ways to improve sensitivity with further filtering and correction.

rkmh is threaded using OpenMP but the code should be considered minimally tuned. Hashing can handle around 250 reads/second and kmer comparison around 40 reads / second.

### Getting help
Please post to the [github](https://github.com/edawson/rkmh.git) for help.
