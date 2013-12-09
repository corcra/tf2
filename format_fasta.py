#!/usr/bin/python

import gzip
import sys


fasta_file = gzip.open(sys.argv[1],'r')
outfile = open(sys.argv[1]+'.oneperline','w')
buffer = ''
for line in fasta_file:
    if '>' in line:
        if len(buffer)!=0:
            outfile.write(buffer.rstrip(' ')+'\n')
        buffer = ''
    else:
        buffer = buffer +' '.join(line.strip().replace('A','1').replace('C','2').replace('G','3').replace('T','4'))+' '
