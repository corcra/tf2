#!/usr/bin/python
# This just takes a (list of) range(s) and makes a fake bedfile with one line per location...

def get_peak_range(pre_line):
    line = pre_line.strip()
    colon = line.find(':')
    hyphen = line.find('-')
    chro = line[:colon]
    start = int(line[(colon+1):hyphen])
    end = int(line[(hyphen+1):])
    return (chro,start,end)

import sys
import math

print 'Expanding peaks!'

if len(sys.argv<3):
    sys.exit('Requires list of peaks and window size!')

peakfile = sys.argv[1]
window_size = float(sys.argv[2])
win_edge = int(math.ceil(window_size/2))
outfile = open(peakfile+'.exp','w')

for line in open(peakfile,'r'):
    [chro,start,end] = get_peak_range(line)
    for i in range(start-win_edge,end+win_edge):
        outline=chro+"\t"+str(i)+"\t"+str(i+1)+"\n"
        outfile.write(outline)
