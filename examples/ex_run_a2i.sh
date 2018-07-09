#!/bin/sh

../scripts/gnumaps.pl --genome ./genomes/ref.fa --pair_1 ./reads/ref-wc-edited_1.txt --pair_2 ./reads/ref-wc-edited_2.txt --read_type rna --per_dist 10 --nt_conv a2i --mpi 1 --num_threads 1 --outdir ./exp/toyRNAe
