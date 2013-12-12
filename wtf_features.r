# Functions pertaining to the pre-processing of the DNase-1 data
args<-commandArgs(TRUE)

# --- Functions --- #
# Remember: signal_data has one line per location, and the data we're assessing is the columns...

# above average signal?
get_feature1 <- function(data,peak_start,peak_end){
    # the first 6 columns are just identifier... who needs identifier, right?
    # pretty sure i'm being paranoid here
#    inrange <- (data[,2] >= peak_start)&&(data[,2]<=peak_end)
#    if (mean(inrange)<1){
#        browser()
#        }
    signal_data <- matrix(sapply(data[,5:ncol(data)],as.numeric),nrow=nrow(data))
    signal_data[signal_data==0]<-NA
    peak_mean <- mean(signal_data[,1],na.rm=TRUE)
    f1<-ifelse(rowMeans(signal_data[,2:ncol(signal_data)],na.rm=TRUE)>peak_mean,2,1)
    f1[is.na(f1)]<-0
    return(f1)
    }

# --- Constants and data --- #

signal<-file("dnase_signal_final.bed.gz",open='r')
feature1<-file("feature1.bed",open="w")
peaklist<-read.table("k562_peak_list",as.is=TRUE)
N_PEAKS<-nrow(peaklist)
for (peak in 1:N_PEAKS){
    hor_str <- peaklist[peak,1]
    if (peak%%10000==0){
        cat("Peak",peak,"\n")
        cat(hor_str,"\n")
    }
    # this is just string formatting
    broken_hor_str <- unlist(strsplit(hor_str,":"))
    peak_range <- unlist(strsplit(broken_hor_str[2],"-"))
    peak_start <- as.integer(peak_range[1])
    peak_end <- as.integer(peak_range[2])
    # this appaers to be the range
    peak_length <- diff(c(peak_start,peak_end))-1
    buff <- matrix(scan(signal,sep="\t",what=character(),nlines=peak_length,quiet=TRUE),nrow=peak_length,byrow=T)
    f1<-get_feature1(buff,peak_start,peak_range)
    chro<-broken_hor_str[1]
    cat(chro,"\t",peak_start,"\t",peak_end,"\t",f1,"\n",file=feature1)
}
