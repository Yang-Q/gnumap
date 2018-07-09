#!/bin/sh

../scripts/gnumaps2.pl --genome ./genomes/ref.fa --pair_1 ./reads/seq1_1.fq --pair_2 ./reads/seq1_2.fq --lib_type unstrand --read_type rna --per_dist 10 --num_amb 10 --num_threads 1 --debug 1 --outdir ./exp/v2_normal
