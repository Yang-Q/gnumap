#!/usr/bin/python

'''
long fa head converter
'''

import os, sys
import subprocess as sp

#============================
def check_if_dos_file(file2):
	print 'checking if %s is in dos format...' % file2
	cmd = 'grep -iP -m1 \'\r$\' %s' % file2
	proc = sp.Popen(cmd, stdout=sp.PIPE, shell=True)
	(msg, err) = proc.communicate()
	msg = msg.strip()
	dosF=False
	if msg:
		dosF=True
		print 'yes, [%s] is a dos file format, it contains [CR+LF]' % file2
	print 'done.'
	return dosF
	
#============================
def replace_crlf2lf_dos2unix(dosfile):
	unixfile = dosfile+'.tmp'
	print 'converting dos file[%s] to unix format...' % dosfile
	cmd = 'tr -d \'\r\' < %s > %s' % (dosfile,unixfile)
	#print cmd #debug
	os.system(cmd)
	os.rename(unixfile,dosfile)
	print 'done'

def longHeadFa2short(fa):
	
	sys.stderr.write('converting short head [%s] ...\n' % fa)
	
	fp = open(fa,'r')

	idx=1
	line = fp.readline()
	map = fa+'.map'
	if line[:5] == '>ref1':
		fp.close()
		if os.path.exists(map):
			return map
		else:
			sys.exit('[error]%s does not exists! prepare original %s file\n' % (map,fa))

	fp.seek(0, 0)
	fa2 = fa+'.tmp'
	fp2 = open(fa2,'w')
	fp2b = open(map,'w')
	
	for i in fp:
		if i[0]=='>':
			fp2b.write('%s\tref%d\n' % (i[1:].rstrip(),idx))
			fp2.write('>ref%d\n' % idx)
			idx+=1
		else:
			fp2.write('%s' % i)
	
	fp.close()
	fp2.close()
	fp2b.close()
	
	#backup the original fa file
	cmd = 'gzip %s' % fa
	os.system(cmd)
	os.rename(fa2,fa)
	sys.stderr.write('done.\n')
	
	return map
	
#============================
def check_if_file_valid(file1):
	validF=False
	if os.path.exists(file1):
		statinfo = os.stat(file1).st_size
		if statinfo>0:
			validF=True
	return validF	
	
#============================
def longFaSam2short(fa,sam):
	
	#check if sam is not empty or exists
	if not check_if_file_valid(sam):
		return 0
	
	sys.stderr.write('converting long head [%s] ...\n' % sam)
	map = fa+'.map'
	if not os.path.exists(map):
		sys.exit('[error]%s does not exists!' % map)
	idx2tran = {}
	fp = open(map,'r')
	for i in fp:
		words=i.split('\t')
		idx2tran[words[1].rstrip()]=words[0]
	fp.close()
	
	fp = open(sam,'r')
	sam2 = sam+'.tmp'
	fp2 = open(sam2,'w')
	for i in fp:
		words=i.split('\t')
		if words[2] in idx2tran:
			words[2] = idx2tran[words[2]]
		fp2.write('%s' % '\t'.join(words))
		
	fp.close()
	fp2.close()
	#sys.exit(1) #debug
	os.rename(sam2,sam)
	sys.stderr.write('done.\n')	
	
	return 1