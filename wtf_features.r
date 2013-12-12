# Functions pertaining to the pre-processing of the DNase-1 data
options(warn=-1)

# --- Functions --- #

# above average signal?
get_feature1 <- function(signal_data,peak_mean){
    col_means <- colMeans(signal_data)
    f1<-ifelse(col_means>peak_mean,2,ifelse(col_means>0,1,0))
    return(f1)
}

get_feature2 <- function(signal_data){
    for (loc in 1:ncol(signal_data)){

        }
    

    }

# --- Constants and data --- #

signal<-file("dnase_signal_final.bed.gz",open='r')
feature1<-file("feature1.bed",open="w")
peaklist<-read.table("k562_peak_list",as.is=TRUE)
N_PEAKS<-nrow(peaklist)
for (peak in 1:N_PEAKS){
    hor_str <- peaklist[peak,1]
    if (peak%%1000==0){
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
    peak_mean <- mean(as.numeric(buff[,5]))
    signal_data<- apply(buff[,5:ncol(buff)],1,as.numeric)
    f1<-get_feature1(signal_data,peak_mean)
    chro<-broken_hor_str[1]
    cat(chro,"\t",peak_start,"\t",peak_end,"\t",f1,"\n",file=feature1)
}
