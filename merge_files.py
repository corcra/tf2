#!/usr/bin/python
# This script will merge as many files as you give it, line by line. (so if you have files A,B,C, the first 3 lines in the file will be A1,B1,C1...)
# WARNING: the files must all have the same number of lines! Otherwise you're in serious trouble.

import gzip
import sys

outfile = open('merged_files','w')

firstfile = gzip.open(sys.argv[1])
otherfiles = [gzip.open(arg) for arg in sys.argv[2:]]

for line in firstfile:
    identifier = line.split()[0:2]
    length = len(line.split())
    outfile.write(line)
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
            outfile.write(otherline)
