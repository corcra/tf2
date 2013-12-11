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

def get_flanks(i,buff,WIN_EDGE):
# problem: we want to create flanks based on GENETIC distance
# best case scenario is each location in the buffer corresponds to each location
# in any other scenario, some locations are MISSING -> the buffer contains a larger region
# so the true flanking site is inside buff[(i-WIN_EDGE):(i+WIN_EDGE+1)]
# we just have find it! -> check each element if it's inside...
    current_loc = int(buff[i].split()[1])
    lower_flank = current_loc-WIN_EDGE
    # note inclusive
    upper_flank = current_loc+WIN_EDGE
    flank_range = range(lower_flank,upper_flank+1)
    flank_list = [el.split()[1]+' '+el.split()[2] for el in buff[(i-WIN_EDGE):(i+WIN_EDGE+1)] if int(el.split()[1]) in flank_range]
    # the output of this script is going to look rather weird
#    flank_list = [el.split()[1]+' '+el.split()[2] for el in buff[(i-WIN_EDGE):(i+WIN_EDGE+1)] if in_range(int(el.split()[1]),lower_flank,upper_flank)]
    flanks = ' '.join(flank_list)
    print flank_range
    print flanks
    return flanks

def process_peak(peak_num,peak_start,peak_end,buff,outfile,new_peaks_file):
    # what we have: one line per location (good!), with CHR LOC VAL PEAK_NUM
    # peak_end is the last location in the peak, not the upper bound
    peak_length = peak_end-peak_start+1
    buff_length = len(buff)
    # what if there isn't enough buffer? NOTE: THIS MODIFIES THE PEAK CALLS, SO WE NEED TO USE THE OUTPUT OF THIS SCRIPT TO DEFINE OUR PEAKS FROM NOW ON... (the new_peaks_file)..
    buff_start = int(buff[0].split()[1])
    buff_end = int(buff[-1].split()[1])
    if (peak_start-buff_start)<WIN_EDGE:
        # move up the start of the peak
        peak_start = buff_start+WIN_EDGE
    if (buff_end-peak_end)<WIN_EDGE:
        # move in the end of the peak
        peak_end = buff_end-WIN_EDGE
    # this will be the total signal in the peak... can rewrite the dnase_peaks file to include this, i guess
    peak_total = 0
    for i in range(buff_length):
        element = buff[i]
        [chr,loc,val] = element.split()
        if in_range(int(loc),peak_start,peak_end):
            flanks = get_flanks(i,buff,WIN_EDGE)
            if peak_num==2:
                print buff[1:10]
                print peak_start, peak_end
                print element.split()
                print i
                sys.exit()
            peak_total = peak_total + float(val)
            loc_line = chr+' '+str(loc)+' '+str(peak_num)+' '+flanks+'\n'
            outfile.write(loc_line)
    new_peaks_file.write(chr+'\t'+str(peak_start)+'\t'+str(peak_end)+'\t'+str(peak_total)+'\n') 
    return

if len(sys.argv)<3:
    sys.exit('Requires .bedGraph and peaklist file!')

# each line of this file contains peak range
peaks_file = open(sys.argv[2],'r')
new_peaks_file = open(sys.argv[2]+'.totals','w')
bedgraph = gzip.open(sys.argv[1],'r')
out_file = open(sys.argv[1]+'.processed','w')

first_peak = peaks_file.readline()
[peak_start,peak_end] = get_peak_range(first_peak)

peak_buff = []

line_num = 0
peak_num = 1
read_peaks = 1
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
            temp_line = chro+' '+str(temp_loc)+' '+str(val/length)
            peak_buff.append(temp_line)
    else:
        if len(peak_buff)>0:
            process_peak(peak_num,peak_start,peak_end,peak_buff,out_file,new_peaks_file)
            peak_buff = []
            peak_num = peak_num + 1
        if end>peak_end:
            # get a new peak
            new_peak = peaks_file.readline()
            read_peaks = read_peaks + 1
            # avoid EOF
            if len(new_peak)>0:
                [peak_start,peak_end] = get_peak_range(new_peak)

print 'read',read_peaks-1,'and wrote',peak_num-1
