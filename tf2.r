# Transcription factor binding site identification! With interactions!

# ---- For my implementation: fix the other factors ---- #
cat('Getting binding status from ChIP-seq data!\n')
#bound_from_chip <- as.matrix(read.table('data/binding_mat',header=T))
bound_from_chip <- as.matrix(read.table('fake_data/fake_binding',header=T))

# ---- Constants! ---- #
FACTORS <- colnames(bound_from_chip)
N_FACTORS <- length(FACTORS)
N_PEAKS <- nrow(bound_from_chip)
#N_PEAKS <- 4
N_ITER <- 5
#N_FEATURES <- 1
N_FEATURES <- 3
EM_THRESHOLD <- 0.5
TAU <- 0.2
N_CORES <- 1

# ---- Load functions! ---- #
source('tf2_functions.r')
#.libPaths(new="/home/sh985/R/x86_64-unknown-linux-gnu-library/3.0")
=======
N_CORES <- 8

# ---- Load functions! ---- #
source('tf2_functions.r')
.libPaths(new="/home/sh985/R/x86_64-unknown-linux-gnu-library/3.0")
library(rqhmm)
library(parallel)

# ---- Load Data! ---- #
cat("Getting data!\n")
#fc <- file('data/processed_data.gz',open='r')
fc <- file('fake_data/fd.gz',open='r')
data <- vector("list",N_PEAKS)
for (peak in 1:N_PEAKS){
    buff <- scan(fc,sep=" ",what=numeric(),nlines=(N_FEATURES+1),quiet=TRUE)
    # columns -> number of locations, rows -> number of emission variables (first one will be DNA sequence)
    data[[peak]] <- matrix(buff,nrow=(N_FEATURES+1),byrow=T)
    }
close(fc)
cat("Data loaded!\n")

# ---- Initialise binding status ---- #
# take it from chip!
binding_status <- bound_from_chip
# initialise our target factor (CTCF) to something random...
binding_status[,"CTCF"]<-rbinom(N_PEAKS,1,0.5)
binding_temp <- binding_status
colnames(binding_status)<-FACTORS
colnames(binding_temp)<-FACTORS

# ---- Initialise parameters! ---- # These are all per-factor! Hence lists!
emission_params <- vector("list",N_FACTORS)
transition_params <- vector("list")
motifs<-vector("list")
# per-peak transition params
for (factor in FACTORS){
    transition_params[[factor]] <- vector("list",N_PEAKS)
    motifs[[factor]]<-get_motif(factor)
}

#ctcf_peaks <- which(bound_from_chip[,"CTCF"]==1)
# for testing: only looping over one TF
TEST_FACTORS<-"CTCF"
delta.binding<-vector("numeric")
# ---- The outer loop: 'sample' over binding configurations ---- #
for (iter in 1:N_ITER){
    cat('Iteration',iter,'\n')
    # ---- First middle loop: iterate over each factor! ---- #
    for (factor in TEST_FACTORS){
        cat('Getting binding status for',factor,'\n')
        # get the motif
        pwm <- motifs[[factor]]
        factor_size <- ncol(pwm)
        N_STATES <- factor_size+2
        # initialise emission parameters for this factor
        if (is.null(emission_params[[factor]])){
            cat("No emission parameters saved for ",factor," - making some up!\n")
            emission_params[[factor]] <- matrix(rep(0.5,N_STATES*N_FEATURES),nrow=N_STATES,ncol=N_FEATURES)
        }
        # initialise a blank hmm with the right size... and fixed variables
        factor_hmm <- build_hmm(N_STATES,N_FEATURES,pwm)
        # calculate coincidence of this TF with the rest (used to get a_BF)
        coincidence <- get_interactions(factor,binding_status)
        # when we finish EM we will record binding predictions for all peaks
        all_peaks_bound <- rep(NA,N_PEAKS)
        # initialise the while loop
        delta.ll <- EM_THRESHOLD*2
        ll.old <- -Inf
        EM.iter <- 0
        ll.all <- vector("numeric")
        decrease <- 0
        # ---- Second middle loop: EM! ---- #
        while(abs(delta.ll)>EM_THRESHOLD){
            EM.iter <- EM.iter + 1
            cat("EM iteration",EM.iter,"\n")
            theta_denom <- rep(0,N_STATES)
            theta_numer <- matrix(rep(0,N_STATES*N_FEATURES),nrow=N_STATES,ncol=N_FEATURES)
            # ll over the peaks
            ll_cumulative <- 0

            # parallel section... this will be lengthy
            # ---- Inner 'loop': iterate over peaks! --- #
            peak_results <- mclapply(1:N_PEAKS,eval_peak,factor_hmm,data,transition_params[[factor]],emission_params[[factor]],N_STATES,N_FEATURES,factor,factor_size,binding_status,coincidence,TAU,mc.cores=N_CORES)

            for (peak in 1:N_PEAKS){
                this_peak <- peak_results[[peak]]

                # theta is additive over peaks, remember
                theta <- this_peak$"theta"
                theta_numer <- theta_numer + theta$"theta_numer"
                theta_denom <- theta_denom + theta$"theta_denom"

                # update transition params per-peak
                transition_params[[factor]][[peak]] <- this_peak$"new_trans_params"
                
                # loglik is additive over peaks
                loglik<-this_peak$"loglik"
                ll_cumulative <- ll_cumulative + loglik
            
                # bound or not? we only care after EM has converged
                byn <- this_peak$"bound"
                if (byn==1){
                    print("we think this peak is bound!")
                    browser()
                }
                all_peaks_bound[[peak]] <- this_peak$"bound"
            }
            # update emission parameters after all peaks - if commented out, we're avoiding that (due to numerical issues)
#            emission_params[[factor]] <- theta_numer/theta_denom

            # check how the likelihood has changed...
            ll <- ll_cumulative
            delta.ll <- ll-ll.old
#            print(ll)
            if(delta.ll<0){
                cat("ERROR: likelihood is decreasing!\n")
                print(delta.ll)
                decrease = decrease + delta.ll
            }
            ll.old <- ll
            ll.all <- c(ll.all,ll)
        }
        cat("EM has converged!\n")
#        plot(ll.all,type='l',xlab='Iteration',ylab='Log-likelihood')
        cat('lhood decreased by',decrease,'in total\n')

        binding_temp[,factor] <- all_peaks_bound
    }
    # for the purpose of somehow gauging if convergence is occurring
    # visualise_binding(binding_status)
    # using the infinity norm here, for the... fun?
    delta.binding <- c(delta.binding,norm((binding_status-binding_temp),"I"))

    binding_status <- binding_temp
}

# ---- After iteration: retrieve predictions ---- #
cm <- get_confusion_matrix(binding_status,bound_from_chip,FACTORS)
columns<-cbind("TF","TP","FP","TN","FN","sens","spec")
write(columns,file="cm.txt",ncolumns=7,sep="\t")
for (factor in FACTORS){
    write(c(factor,unlist(cm[[factor]])),file="cm.txt",append=TRUE,sep="\t",ncolumns=7)
}
save.image('tf2.RData')
