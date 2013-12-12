# --- these are the input files --- #
DNASE_PEAKS=wgEncodeOpenChromDnaseK562PkV2.narrowPeak
DNASE_SIGNAL=wgEncodeUwDgfK562Sig.bedGraph.gz
WIN_SIZE=50

echo "Getting peaks!"
awk '{ print $1":"$2"-"$3+1 }' $DNASE_PEAKS> k562_peak_list

# expand the peaks! ... creates k562_peak_list.exp
python expand_peak.py k562_peak_list $WIN_SIZE
echo "Restricting DNAse signal to be in/near peaks!"
zcat $DNASE_SIGNAL | bedmap --echo --skip-unmapped - k562_peak_list.exp > signal_near_peaks.bed
gzip signal_near_peaks.bed
# expand thie bedgraph! ... this generates signal_near_peaks.bed.exp.gz
python expand_bedgraph.py signal_near_peaks.bed.gz
echo "Assigning scores to locations in peaks!"
zcat signal_near_peaks.bed.exp.gz | bedmap --exact --echo --echo-map-score --delim '\t' k562_peak_list.exp - > prescored.bed
# fill in zeroes!
awk '{ if (NF==3) print $1"\t"$2"\t"$3"\tNA\t"0; else print $1"\t"$2"\t"$3"\tNA\t"$4 }' prescored.bed > scored_in_peaks.bed
# get the flanking regions for each location (overlapping with self!)
echo "Getting flanking signals!"
bedmap --range 25 --echo --echo-map-score scored_in_peaks.bed | sed 's/[;|\|]/\t/g' > signal_with_flanks.bed
# now re-restrict this to the original peak regions
echo "Re-restricting locations (with flank info) to be inside peaks!"
bedmap --echo --skip-unmapped signal_with_flanks.bed $DNASE_PEAKS  > dnase_signal_final.bed
# ready to go to R!
gzip dnase_signal_final.bed
# clean up!
echo "Cleaning up!"
rm -v k562_peak_list.exp
rm -v signal_near_peaks.bed.gz
rm -v signal_near_peaks.bed.exp.gz
rm -v prescored.bed
rm -v scored_in_peaks.bed
rm -v signal_with_flanks.bed
# get the features!
echo "Extracting features with R!"
R --file=extract_features.r --args dnase_signal_final.gz

# next part... sequence
echo "Getting sequence!"
awk '{ print $1":"$2"-"$3 }' $DNASE_PEAKS > peaklist_for_seq
twoBitToFa /gbdb/hg19/hg19.2bit peakseq.fa -seqList=peaklist_for_seq -noMask
gzip peakseq.fa
# format the sequence! ... this generates peakseq.fa.gz.oneperline
python format_fasta.py peakseq.fa.gz
# overlap the sequence with the peaks...
echo "Overlapping sequence!"
bedmap --echo --skip-unmapped peakseq.fa.gz.oneperline $DNASE_PEAKS > seq.bed
gzip seq.bed
# clean up
rm -v peakseq.fa.gz.oneperline

# merge! this produces processed_data.gz which can be fed into the program...
echo "Merging sequence and DNase features!"
python merge_files.py seq.bed.gz feature1.bed.gz feature2.bed.gz ... featureN.bed.gz
# final clean up
# rm -v feature1.bed.gz
# ...
# rm -v featureN.bed.gz
