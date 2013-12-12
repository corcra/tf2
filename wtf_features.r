# Functions pertaining to the pre-processing of the DNase-1 data
args<-commandArgs(TRUE)

# --- Functions --- #
# Remember: signal_data has one line per location, and the data we're assessing is the columns...

# above average signal?
get_feature1 <- function(data,peak_start,peak_end){
    # the first 6 columns are just identifier... who needs identifier, right?
    inrange <- (data[,2] >= peak_start)&&(data[,2]<=peak_end)
    if (mean(inrange)<1){
        browser()
        }
    signal_data <- data[,6:ncol(data)]
    f1<-1
#    f1<-ifelse(rowMeans(signal_data)>mean,2,1)
    return(f1)
    }

# --- Constants and data --- #

signal<-file(args[1],open='r')
feature1<-file("feature1.bed",open="w")
peaklist<-read.table(args[2])
N_PEAKS<-nrow(peaklist)
for (peak in 1:N_PEAKS){
    hor_str <- peaklist[peak]
    cat(hor_str,"\n")
    # this is just string formatting
    peak_range <- unlist(lapply(strsplit(unlist(strsplit(hor_str,":"))[2],"-"),as.integer))
    peak_start <- peak_range[1]
    peak_end <- peak_range[2]
    peak_length <- diff(peak_range)+1
    buff <- matrix(scan(signal,sep="\t",what=character(),nlines=peak_length,quiet=TRUE),nrow=peak_length,byrow=T)
    f1<-get_feature1(buff,peak_start,peak_range)
    cat(chro,"\t",peak_start,"\t",peak_end,"\t",f1,"\n",file=feature1)
    }

