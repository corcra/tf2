#!/usr/bin/python
# All this does is rewrite a fasta file so that the whole sequence is on one line, there's no identifier (perhaps dubious, might change this) and recoding ACGT -> 1234 (for use in the HMM later!)

import gzip
import sys

def get_peak_range(pre_line):
    line = pre_line.strip()
    colon = line.find(':')
    hyphen = line.find('-')
    chro = line[1:colon]
    start = int(line[(colon+1):hyphen])
    end = int(line[(hyphen+1):])
    return (chro,start,end)

fasta_file = gzip.open(sys.argv[1],'r')
outfile = open(sys.argv[1]+'.oneperline','w')
buff = ''
ignored_missing = 0
for line in fasta_file:
    if '>' in line:
        if len(buff)!=0:
            outfile.write(buff.rstrip(' ')+'\n')
        [chro,start,end] = get_peak_range(line)
        buff = chro+'\t'+str(start)+'\t'+str(end)+'\t'
    else:
        buff = buff +' '.join(line.strip().replace('A','1').replace('C','2').replace('G','3').replace('T','4')).replace('N','0')+' '

