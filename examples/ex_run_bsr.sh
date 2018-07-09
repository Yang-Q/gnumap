#!/bin/sh

../scripts/gnumaps.pl --genome ./genomes/ref1k.fa --pair_1 ./reads/ref-wc-methy_1.txt --pair_2 ./reads/ref-wc-methy_2.txt --read_type dna --per_dist 10 --nt_conv bs --num_threads 1 --pileup 1 --outdir ./exp/toyBSR
