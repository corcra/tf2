    # get the peaks!
zcat k562_dnase/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz | awk '{ print $1":"$2"-"$3+1 }' > k562_dnase_peaks
    # get the sequence!
twoBitToFa /gbdb/hg19/hg19.2bit peak_seq.fa -seqList=k562_dnase_peaks -noMask
gzip peak_seq.fa
    # format the sequence!
python format_fasta.py peak_seq.fa.gz
    # overlap with called peaks to get peak VALUES (in bedgraph format)!
zcat ~/k562_dnase/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz | bedmap --range 100 --echo --skip-unmapped wgEncodeOpenChromDnaseK562SigV2.bedGraph - > signal_in_peaks.bed
gzip signal_in_peaks.bed
    # process this into R-ready format! warning: this takes ages. second warning: file locations... check them!
python format_bedgraph.py signal_in_peaks.bed.gz
gzip processed_signal_in_peaks.bed

    # SPECIFIC TO MY DATA: the above bedgraph: signal_in_peaks appears to be a subset of the sequence data! ... except for a single region on the X chromosome? (not sure what's going on here!)
    # THEREFORE: we shall trim down the sequence data so that it only contains those lines also in the signal data... (at this point, individual lines refer to peaks)... format_bedgraph produces a list of recorded peaks, so we use this
bedmap --echo --skip-unmapped peak_seq.fa.gz.oneperline k562_dnase_peaks_totals > sequence_by_peaks.bed
gzip sequence_by_peaks.bed
    # clean up
rm peak_seq.fa.gz.oneperline

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
