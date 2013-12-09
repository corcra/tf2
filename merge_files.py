#!/usr/bin/python
# This script will merge as many files as you give it, line by line. (so if you have files A,B,C, the first 3 lines in the file will be A1,B1,C1...)
# WARNING: the files must all have the same number of lines! Otherwise you're in serious trouble.

import gzip
import sys

outfile = open('merged_files','w')

firstfile = gzip.open(sys.argv[1])
otherfiles = [gzip.open(arg) for arg in sys.argv[2:]]

for line in firstfile:
    outfile.write(line)
    for otherfile in otherfiles:
        otherline = otherfile.readline()
        outfile.write(otherline)
