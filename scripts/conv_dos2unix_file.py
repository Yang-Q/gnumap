#!/usr/bin/python

'''
file format conversion
'''

import os, sys
import gnumaps_lib

fmt2check=sys.argv[1]
if gnumaps_lib.check_if_dos_file(fmt2check):
	gnumaps_lib.replace_crlf2lf_dos2unix(fmt2check)