#!/usr/bin/python
# This script will merge as many files as you give it, line by line. (so if you have files A,B,C, the first 3 lines in the file will be A1,B1,C1...)
# WARNING: the files must all have the same number of lines! Otherwise you're in serious trouble.

import gzip
import sys

outfile = gzip.open('processed_data.gz','w')

firstfile = gzip.open(sys.argv[1])
otherfiles = [gzip.open(arg) for arg in sys.argv[2:]]

line_num = 0
for line in firstfile:
    if line_num%10000==0:
        print line_num
    line_num = line_num+1
    identifier = line.split()[0:2]
    length = len(line.split())
    outfile.write(' '.join(line.split()[3:])+'\n')
    for otherfile in otherfiles:
        otherline = otherfile.readline()
        other_ident = otherline.split()[0:2]
        other_length = len(otherline.split())
        if not other_ident == identifier:
            print 'ERROR! DIFFERING IDENTIFIERS?'
            print identifier, other_ident
            sys.exit()
        elif not other_length == length:
            print 'ERROR! DIFFERING LENGTHS?'
            print length, other_length
            print line.split()
            print otherline.split()
            sys.exit()
        else:
            outfile.write(' '.join(otherline.split()[3:])+'\n')
