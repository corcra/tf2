# --- these are the input files --- #
DNASE_PEAKS=wgEncodeOpenChromDnaseK562PkV2.narrowPeak
DNASE_SIGNAL=wgEncodeUwDgfK562Sig.bedGraph.gz
WIN_SIZE=50

echo "Getting peaks!"
awk '{ print $1":"$2"-"$3+1 }' $DNASE_PEAKS> k562_peak_list

# expand the peaks! ... this creates $PEAK_LIST.exp (by the way, this produces 'expanded' peaks, including the flanking regions, so they need to be re-restricted later!
python expand_peak.py k562_peak_list $WIN_SIZE
echo "Restricting DNAse signal to be in/near peaks!"
zcat $DNASE_SIGNAL | bedmap --echo --skip-unmapped - k562_peak_list.exp > signal_near_peaks.bed
gzip signal_near_peaks.bed
# expand this bedgraph (one score per location...) ... this generates $SIGNAL_IN_PEAKS.exp.gz
python expand_bedgraph.py signal_near_peaks.bed.gz
echo "Assigning scores to locations in peaks!"
zcat signal_near_peaks.bed.exp.gz | bedmap --exact --echo --echo-map-score --delim '\t' k562_peak_list.exp - > prescored.bed
# fill in zeroes!
awk '{ if (NF==3) print $1"\t"$2"\t"$3"\tNA\t"0; else print $1"\t"$2"\t"$3"\tNA\t"$4 }' prescored.bed > scored_in_peaks.bed
# now get the flanking regions for each location in each peak (the sed stuff is for getting rid of semicolons) (overlapping with self!)
echo "Getting flanking signals!"
bedmap --range 25 --echo --echo-map-score scored_in_peaks.bed | sed 's/[;|\|]/\t/g' > signal_with_flanks.bed
# now re-restrict this to the original peak regions
echo "Re-restricting locations (with flank info) to be inside peaks!"
bedmap --echo --skip-unmapped signal_with_flanks.bed $DNASE_PEAKS  > dnase_signal_final.bed
# ready to go to R!
gzip dnase_signal_final.bed
# clean up!
echo "Cleaning up!"
rm -v prescored.bed
rm -v scored_in_peaks.bed
rm -v signal_with_flanks.bed
rm -v signal_near_peaks.bed.exp.gz
rm -v signal_near_peaks.bed.gz
#
exit
# get the features!
echo "Extracting featuers with R!"
R --file=wtf_features.r --args dnase_signal_final.gz

# get the sequence!
awk '{ print $1":"$2"-"$3 }' $DNASE_PEAKS > peaklist_for_seq
twoBitToFa /gbdb/hg19/hg19.2bit peakseq.fa -seqList=peaklist_for_seq -noMask
gzip peakseq.fa
# format the sequence! ... this creates a file called $PEAK_SEQ.gz.oneperline and also $PEAK_LIST.totals
python format_fasta.py peakseq.fa.gz

# SPECIFIC TO MY DATA: the above bedgraph: signal_in_peaks appears to be a subset of the sequence data! ... except for a single region on the X chromosome? (not sure what's going on here!)
# THEREFORE: we shall trim down the sequence data so that it only contains those lines also in the signal data... (at this point, individual lines refer to peaks)... format_bedgraph produces a list of recorded peaks, so we use this
bedmap --echo --skip-unmapped $PEAK_SEQ.gz.oneperline $DNASE_PEAKS > seq.bed
gzip seq.bed
# clean up
rm -v peakseq.fa.gz.oneperline

# by now we should have
# signal_in_peaks and sequence_by_peaks, both gzipped 'bed' files...
# sequence_by_peaks is ONE LINE PER PEAK, ready to be merged
# signal_in_peaks is ONE LINE PER LOCATION, requiring additional formatting
# Also: k562_dnase_peaks_totals which is ONE LINE PER PEAK
# after formatting signal_in_peaks (into separate files per feature - one line per peak) we need to overlap with sequence_by_peaks on each of the resultant feature files to make sure they're 'aligned', as it were
# then the merging can occur!
    
# Re the formatting signal_in_peaks: this requires...
# Read the output of format_bedgraph... this has one line per location, and a list of the DNase values around that location!
# These need to be processed into features! (probably use R to do this)
# For each feature:
# Output a new file with one line PER PEAK: columns are feature values for the locations in the peak.
# Use merge_files.py on the seq data (output from format_fasta.py) and ALL the feature files!
# There are probably better ways to do this (there are definitely better ways to do this) but this can do for now...

    # doing the above ... this will produce merged_files.gz, which can then be used as input to the algorithm!
python merge_files.py peak_seq_fa.gz.oneperline.gz feature1.gz feature2.gz ... featureN.gz



# --- working files --- #
SIGNAL_IN_PEAKS=signal_in_peaks.bed
PEAK_LIST=k562_dnase_peaks
PEAK_SEQ=peak_seq.fa
PEAK_SEQ_WITHSIGNAL=seq_in_peaks.bed
CHIP_PEAKS=peakSeq.Haibk562Ctcf.bed.gz


