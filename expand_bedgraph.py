#!/usr/bin/python

# This will take a bedgraph with CHR START END SCORE and expand all elements (eg if END-START>1)...

import sys
import gzip

if not '.bed.gz' in sys.argv[1]:
    sys.exit('Requires gzipped bed file!')

print 'Expanding bedgraph!'

outfile = gzip.open(sys.argv[1].replace('.bed.gz','.bed.exp.gz'),'w')

linenum = 0
for line in gzip.open(sys.argv[1],'r'):
    if linenum%50000==0:
        print linenum
    linenum=linenum+1
    [chro,start,end,score] = line.split()
    length = int(end)-int(start)
    for n in range(length):
        newloc = int(start)+n
        sc = float(score)/length
        newline = chro+'\t'+str(newloc)+'\t'+str(newloc+1)+'\tNA\t'+str(sc)+'\n'
        outfile.write(newline)
