#!/usr/bin/python
# This will take a bedgraph and expand it into a weird sort of bed file, where we have
# CHR LOC DNASE_1 ... DNASE_N
# N is the size of the window around LOC, and DNASE_i are the DNASE values corresponding the window locations
# The output of this script shall be the input to the feature extraction part of the larger code
import gzip
import sys

# arbitrary!
WIN_SIZE = 25
WIN_EDGE = WIN_SIZE/2

def overlap(A1,A2,B1,B2):
    if A1 <= B1 and B1 <= A2:
        return True
    elif B1 <= A1 and A1 <= B2:
        return True
    else:
        return False

def in_range(x,A,B):
    if x < A or x > B:
        return False
    else:
        return True

def get_peak_range(pre_line):
    line = pre_line.strip()
    colon = line.find(':')
    hyphen = line.find('-')
    start = int(line[(colon+1):hyphen])
    end = int(line[(hyphen+1):])
    return (start,end)

def get_flanks(i,buffer,WIN_EDGE):
    flank_list = [element.split()[2] for element in buffer[(i-WIN_EDGE):(i+WIN_EDGE+1)]]
    flanks = ' '.join(flank_list)
    return flanks

def process_peak(peak_start,peak_end,buffer,outfile,new_peaks_file):
    # what we have: one line per location (good!), with CHR LOC VAL
    # peak_end is the last location in the peak, not the upper bound
    peak_length = peak_end-peak_start+1
    buffer_length = len(buffer)
    # this will be the total signal in the peak... can rewrite the dnase_peaks file to include this, i guess
    peak_total = 0
    for i in range(buffer_length):
        element = buffer[i]
        [chr,loc,val] = element.split()
        if in_range(int(loc),peak_start,peak_end):
            peak_total = peak_total + float(val)
            flanks = get_flanks(i,buffer,WIN_EDGE)
            loc_line = chr+' '+str(loc)+' '+flanks+'\n'
            outfile.write(loc_line)
    new_peaks_file.write(chr+'\t'+str(peak_start)+'\t'+str(peak_end)+'\t'+str(peak_total)+'\n') 

# each line of this file contains peak range
peaks_file = open('../k562_dnase_peaks','r')
new_peaks_file = open('../k562_dnase_peaks_totals','w')
bedgraph = gzip.open(sys.argv[1],'r')
out_file = open('../dnase_in_windows','w')

first_peak = peaks_file.readline()
[peak_start,peak_end] = get_peak_range(first_peak)

peak_buffer = []

line_num = 0
for line in bedgraph:
    # file's pretty big
    if line_num%100000==0:
        print line_num
    line_num = line_num+1

    chro = line.split()[0]
    start = int(line.split()[1])
    end = int(line.split()[2])
    if overlap(start,end,(peak_start-WIN_EDGE),(peak_end+WIN_EDGE)):
        length = end-start
        val = float(line.split()[3])
        # creating a buffer of all the locs (expanded - one loc per line) in the peak
        for n in range(length):
            temp_loc = start+n
            temp_line = chro+' '+str(temp_loc)+' '+str(val)
            peak_buffer.append(temp_line)
    else:
        if len(peak_buffer)>0:
            process_peak(peak_start,peak_end,peak_buffer,out_file,new_peaks_file)
            peak_buffer = []
        if end>peak_end:
            # get a new peak
            new_peak = peaks_file.readline()
            # avoid EOF
            if len(new_peak)>0:
                [peak_start,peak_end] = get_peak_range(new_peak)

