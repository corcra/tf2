# Functions pertaining to the pre-processing of the DNase-1 data
options(warn=-1)

# --- Functions --- #

# above average signal?
get_feature1 <- function(signal_data,peak_mean){
    col_means <- colMeans(signal_data)
    f1<-ifelse(col_means>peak_mean,2,1)
    return(f1)
}

get_f2_and_f3 <- function(signal_data){
    win_size <- nrow(signal_data)
    # centred
    f2<-vector("integer")
    f3<-vector("integer")
    x <- seq(win_size)-ceiling(win_size/2)
    for (loc in 1:ncol(signal_data)){
        linear <- lm(signal_data[,loc] ~ x)
        inc <- (linear$coefficients["x"]>1/win_size)*1
        quad <- lm(signal_data[,loc] ~ x^2)
        quadratic <- abs(quad$coefficients["x"]>0.3)*1
        f2<-c(f2,inc)
        f3<-c(f3,quadratic)
    }
    return(list("f2"=f2,"f3"=f3))
}

# --- Constants and data --- #

signal<-file("dnase_signal_final.bed.gz",open='r')
#feature1<-file("feature1.bed",open="w")
feature2<-file("feature2.bed",open="w")
feature3<-file("feature3.bed",open="w")
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
    # signal data is now actually columns for locations and rows for the values...
    signal_data<- apply(buff[,5:ncol(buff)],1,as.numeric)
#    f1<-get_feature1(signal_data,peak_mean)
    f2_and_f3<-get_f2_and_f3(signal_data)
    chro<-broken_hor_str[1]
    f2<-f2_and_f3$"f2"
    f3<-f2_and_f3$"f3"
#    cat(chro,"\t",peak_start,"\t",peak_end,"\t",f1,"\n",file=feature1)
    cat(chro,"\t",peak_start,"\t",peak_end,"\t",f2,"\n",file=feature2)
    cat(chro,"\t",peak_start,"\t",peak_end,"\t",f3,"\n",file=feature3)
}
