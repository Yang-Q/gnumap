#!/usr/bin/python
'''
postprocessing of RNA-Seq based on REDItools
Now, I simply copy and pasted the example described in the reditools homepage. Later, I will be able to conver to it to OO
'''
import sys, os, re, random
import subprocess as sp
import numpy as np

#class Reditools:
	#dnaRnaBin=''
	#rnaBam
	#dnaBam
	#ref
	#threads
	#outPrefix
	#intFolder
	#minCvg
	#fqOffset
	#minQual
	#minMap
	#minHpoly
	#strand
	
def reditools_post(redit_parser):
	a=redit_parser.parse_args()
	
	method='DnaRna'
	if args.filter_method != 'dna_rna':
		system.exit(1)
		method='denovo'
	
	resultPrefix='complete'
	#cmdOpt='-c 10,1 -Q 33,64 -q 25,25 -m 20,20 -s 2 -S -g 1 -u -a 6-0 -v 2 -n0.0 -N0.0 -V'
	cmdOpt='-F %s -c 10,2 -m 20,20 -s 2 -g 1 -d -e -p -u -a 6-0 -v 1 -N 0.0 -n 0.0 -t 23' % resultPrefix #used in Picardi et. al.
	cmd='%s/REDItoolDnaRna.py -i %s -j %s -f %s -o %s %s' % (a.reditD, a.rna_aln, a.dna_aln, a.ref, a.out_prefix, cmdOpt)
	print cmd
	os.system(cmd)
	outTab='%s/%s_%s/outTable_%s' % (a.out_prefix,method,resultPrefix,resultPrefix)
	
	optSelect='-d 12 -c 2 -C 10 -v 2 -V 0 -f 0.1 -F 1.0 -e -u'
	candFn='%s/%s_%s/candidates.txt' % (a.out_prefix,method,resultPrefix)
	cmd='%s/selectPositions.py -i %s %s -o %s' % (a.reditD, outFn, optSelect, candFn)
	print cmd
	os.system(cmd)
	
	#annotate ALU
	candMaskFn='%s/%s_%s/candidates.rmsk.txt' % (a.out_prefix,method,resultPrefix)
	cmd='%s/AnnotateTable.py -a %s -i %s -u -c 1,2,3 -n RepMask -o %s' % (a.reditD, a.rmsk, candFn, candMaskFn)
	print cmd
	os.system(cmd)
	
	#annotate RefSeq
	candAnnFn='%s/%s_%s/candidates.rmsk.ann.txt' % (a.out_prefix,method,resultPrefix)
	cmd='%s/AnnotateTable.py -a %s -i %s -u -c 1,2 -n RefSeq -o %s' % (a.reditD, a.refgen, candAnnFn)
	print cmd
	os.system(cmd)

#===================================================
#main()
parser = argparse.ArgumentParser(description="rna-edit read simulator [cjhong@bu.edu]")
parser.add_argument('-method', action='store', dest='filter_method', required=True, default='dna_rna', help='[dna_rna] or denovo; if you have both dna and rna alignment file, use dna_rna')
parser.add_argument('-dna', action='store', dest='dna_aln', required=False, help='-')
parser.add_argument('-rna', action='store', dest='rna_aln', required=True, help='-')
parser.add_argument('-ref', action='store', dest='ref', required=True, help='-')
parser.add_argument('-out', action='store', dest='out_prefix', required=True, help='-')
parser.add_argument('-R', action='store', dest='reditD', required=True, help='specify REDItools installation directory')

parser.add_argument('-refgen', action='store', dest='refgen', required=True, help='specify refGen file.')
parser.add_argument('-rmsk', action='store', dest='rmsk', required=True, help='specify rmask annotation file')

parser.add_argument('--reuse', action='store_const', dest='reuse', required=False, const=True, default=False, help='set this flag if you want to run in debugging mode')

#==========
#check if this is DnaRna method or denovo method

#check input option requirement

reditools_post(parser)
