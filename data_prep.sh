DNASE_PEAKS=wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz
DNASE_SIGNAL=wgEncodeOpenChromDnaseK562SigV2.bedGraph
SIGNAL_IN_PEAKS=signal_in_peaks.bed
PEAK_LIST=k562_dnase_peaks
PEAK_SEQ=peak_seq.fa
PEAK_SEQ_WITHSIGNAL=seq_in_peaks.bed
CHIP_PEAKS=peakSeq.Haibk562Ctcf.bed.gz


# get the peaks!
zcat $DNASE_PEAKS | awk '{ print $1":"$2"-"$3+1 }' > $PEAK_LIST
# get the sequence!
twoBitToFa /gbdb/hg19/hg19.2bit $PEAK_SEQ -seqList=$PEAK_LIST -noMask
gzip $PEAK_SEQ
# format the sequence! ... this creates a file called $PEAK_SEQ.gz.oneperline and also $PEAK_LIST.totals
python format_fasta.py $PEAK_SEQ.gz
# overlap with called peaks to get peak VALUES (in bedgraph format)!
zcat $DNASE_PEAKS | bedmap --range 100 --echo --skip-unmapped  $DNASE_SIGNAL - > $SIGNAL_IN_PEAKS
gzip $SIGNAL_IN_PEAKS
# process this into R-ready format! warning: this takes ages. second warning: file locations... check them! ... this creates a file called $SIGNAL_IN_PEAKS.processed ... also $PEAK_LIST.totals ... both of these need to be read into R!
python format_bedgraph.py $SIGNAL_IN_PEAKS.gz $PEAK_LIST
# this file is a weird format, one line per location... the result of the R processing will be one line per peak!
gzip $SIGNAL_IN_PEAKS.gz.processed

# SPECIFIC TO MY DATA: the above bedgraph: signal_in_peaks appears to be a subset of the sequence data! ... except for a single region on the X chromosome? (not sure what's going on here!)
# THEREFORE: we shall trim down the sequence data so that it only contains those lines also in the signal data... (at this point, individual lines refer to peaks)... format_bedgraph produces a list of recorded peaks, so we use this
bedmap --echo --skip-unmapped $PEAK_SEQ.gz.oneperline $PEAK_LIST.totals > $PEAK_SEQ_WITHSIGNAL
gzip $PEAK_SEQ_WITHSIGNAL
# clean up
rm $PEAK_SEQ.gz.oneperline

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



# -- Now for ChIP-seq! -- #
#bigBedToBed raw_chip bedchip...
# Need to look for ChIP-seq peaks in the DNase peaks! the peak list is $PEAK_LIST.totals ... alternately $PEAK_SEQ_WITHSIGNAL.gz should contain this information (these should agree on peaks... need to double check this! maintain consistency!)
zcat $CHIP_PEAKS | bedmap --echo --indicator --delim '\t' $PEAK_LIST.totals - > $FACTOR_BOUND_PEAKS
