#!/usr/bin/python

'''
long fa head converter
'''

import os, sys
import gnumaps_lib

mode=int(sys.argv[1])
fa=sys.argv[2]
sam=sys.argv[3]

if mode == 1: # long head FA --> short head FA
	map = gnumaps_lib.longHeadFa2short(fa)
else: #short head SAM --> long head SAM
	outcode = gnumaps_lib.longFaSam2short(fa,sam)
