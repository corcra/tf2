# extract_features.r - part of TF2
# Purpose: Take appropriately formatted (what is that?) DNase-1 signal and 'extract features' (binary). 0 means 'missing data', 1 'no', 2 'yes'.
# Warning: parallelising this means the peaks might (requires testing) get all messed up, order-wise, so you need to sort-bed on this afterwards! Hopefully this doesn't negate the speed gain from parallelisation.
# Possible solution: mclapply returns the value for each peak, then write this all to the files...


options(warn=-1)
# missingness threshold: what frac of 0s in the signal counts as missing?
THRE<-0.2

# above average signal?
get_feature1 <- function(signal_data,peak_mean){
    col_means <- colMeans(signal_data)
    f1<-ifelse(col_means>peak_mean,2,1)
    f1[is.na(f1)]<-0
    return(f1)
}

get_f2_and_f3 <- function(signal_data){
    win_size <- nrow(signal_data)
    f2<-vector("integer")
    f3<-vector("integer")
    # centred
    x <- seq(win_size)-ceiling(win_size/2)
    for (loc in 1:ncol(signal_data)){
        if (mean(signal_data[,loc]==0)>THRE){
            # missing data
            increasing<-0
            quadratic<-0
        } else{
            linear <- lm(signal_data[,loc] ~ x)
            increasing <- ifelse(linear$coefficients["x"]>1/win_size,2,1)
            quad <- lm(signal_data[,loc] ~ x^2)
            # note presence of another arbitrary threshold
            quadratic <- ifelse(abs(quad$coefficients["x"])>0.3,2,1)
        }
        f2<-c(f2,increasing)
        f3<-c(f3,quadratic)
    }
    return(list("f2"=f2,"f3"=f3))
}

get_peak_features <- function(peak,peaklist,signal,feature1,feature2,feature3){
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
    peak_length <- diff(c(peak_start,peak_end))-1
    buff <- matrix(scan(signal,sep="\t",what=character(),nlines=peak_length,quiet=TRUE),nrow=peak_length,byrow=T)
    peak_mean <- mean(as.numeric(buff[,5]))
    # signal data is now actually columns for locations and rows for the values...
    signal_data<- apply(buff[,5:ncol(buff)],1,as.numeric)
    f1<-get_feature1(signal_data,peak_mean)
    f2_and_f3<-get_f2_and_f3(signal_data)
    chro<-broken_hor_str[1]
    f2<-f2_and_f3$"f2"
    f3<-f2_and_f3$"f3"
    cat(chro,"\t",peak_start,"\t",peak_end,"\t",f1,"\n",file=feature1)
    cat(chro,"\t",peak_start,"\t",peak_end,"\t",f2,"\n",file=feature2)
    cat(chro,"\t",peak_start,"\t",peak_end,"\t",f3,"\n",file=feature3)
}

# --- Main bit! --- #
#
signal<-file("dnase_signal_final.bed.gz",open="r")
feature1<-file("feature1.bed",open="w")
feature2<-file("feature2.bed",open="w")
feature3<-file("feature3.bed",open="w")
peaklist<-read.table("k562_peak_list",as.is=TRUE)
N_PEAKS<-nrow(peaklist)
# i can already feel the race conditions
mclapply(1:N_PEAKS,get_peak_features,peaklist,signal,feature1,feature2,feature3)
