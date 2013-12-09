#!/usr/bin/python
# All this does is rewrite a fasta file so that the whole sequence is on one line, there's no identifier (perhaps dubious, might change this) and recoding ACGT -> 1234 (for use in the HMM later!)

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
