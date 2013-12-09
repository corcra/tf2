# Functions pertaining to the pre-processing of the DNase-1 data
# Input data is the output of format_bedgraph.py... (currently 'pre_process_example') ... also take the list of peaks!

# --- Functions --- #
# Remember: signal_data has one line per location, and the data we're assessing is the columns...

# above average signal?
get_feature1 <- function(signal_data,fc,mean){
    f1<-ifelse(rowMeans(signal_data)>mean,2,1)
    return(f1)
    }

#get_feature2 <- function(signal_data,fc){
#    return(c(2,2,2,1))   
#    }

#get_feature3 <- function(signal_data,fc){
#    return(c(1,2,1,2))   
#    }
#
#
# --- Constants and data --- #

signal<-read.table("processed_signal_in_peaks.bed.gz",as.is=TRUE)
#signal<-read.table("pre_process_example",as.is=TRUE)
peaklist<-read.table("k562_dnase_peaks_totals")
#N_PEAKS <- nrow(peaklist)

# make these manually since I have to specify the features anyway...
# REMEMBER! first 3 fields are identifier, one peak per line, one column per location in the peak!
feature1<-file("feature1.bed",open="w")
#feature2<-file("feature2.bed",open="w")
#feature3<-file("feature3.bed",open="w")

# split the signal into peaks, based on their peak numbers... (third column)
split_signal <- split(signal,list(signal$V3))

# --- Main loopy bit --- #
for (peak in names(split_signal)){
    cat("Peak",peak,"\n")
    peak_signal <- split_signal[[peak]]
    peak_length <- nrow(peak_signal)
    peak_total <- peaklist[peak,4]
    peak_mean <- peak_total/peak_length

    peak_start <- range(peak_signal[,2])[1]
    peak_end <- range(peak_signal[,2])[2]
    chro <- peak_signal[1,1]
    print(chro)

    f1<-get_feature1(peak_signal[,4:peak_length],peak_mean)
#    f2<-get_feature2(peak_signal[,4:peak_length])
#    f3<-get_feature3(peak_signal[,4:peak_length])
    cat(chro,"\t",peak_start,"\t",peak_end,"\t",f1,"\n",file=feature1)
#    cat(chro,"\t",peak_start,"\t",peak_end,"\t",f2,"\n",file=feature2)
#    cat(chro,"\t",peak_start,"\t",peak_end,"\t",f3,"\n",file=feature3)
    }
