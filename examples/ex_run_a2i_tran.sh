#!/bin/sh -l

export WRKD=/unprotected/projects/johnsonlab/exp/cj/gnumaps_tranc_debug/reditDebug

../scripts/gnumaps_pblat.pl --genome $WRKD/genome.fa --transc $WRKD/ce10.transcripts.fasta --pair_1 $WRKD/8374X1b_1.txt --pair_2 $WRKD/8374X1b_2.txt --phred_offset 64 --read_type rna --per_dist 10 --nt_conv a2i --num_threads 8 --outdir $WRKD/out3 --debug 1 
