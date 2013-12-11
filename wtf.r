# Transcription factor binding site identification! Including interactions!

# ---- For my implementation: fix the other factors ---- #
cat('Getting binding status from ChIP-seq data!\n')
bound_from_chip <- as.matrix(read.table('chip_binding_mat',header=T))

# ---- Constants! ---- #
FACTORS <- colnames(bound_from_chip)
N_FACTORS <- length(FACTORS)
N_PEAKS <- nrow(bound_from_chip)
N_ITER <- 5
N_FEATURES <- 1
EM_THRESHOLD <- 0.5
TAU <- 0.1

# ---- Load functions! ---- #
source('wtf_fns_ok.r')
source('wtf_fns.r')
#source('wtf_features.r')
library(rqhmm)

# ---- Process the DNase-1 data! ---- #
# thinking I'll preprocess this...
#DNase_emissions <- get_features(DNase_data)

# ---- Load Data! ---- #
# right now this data is real sequence data but made-up features
fc <- file('processed_data_nomissing.gz',open='r')
data <- vector("list",N_PEAKS)
for (peak in 1:N_PEAKS){
    buff <- scan(fc,sep=" ",what=numeric(),nlines=(N_FEATURES+1),quiet=TRUE)
    # columns -> number of locations, rows -> number of emission variables (first one will be DNA)
    data[[peak]] <- matrix(buff,nrow=(N_FEATURES+1),byrow=T)
    len<-ncol(data[[peak]])
    #data[[peak]][2,] <- rbinom(len,1,0.5)+1
    data[[peak]][2,] <- seq(len)%%2+1
    }
close(fc)
cat("Data loaded!\n")

# ---- Initialise binding status ---- #
# take it from chip!
binding_status <- bound_from_chip
binding_temp <- binding_status
# initialise our target factor (CTCF) to something random...
binding_status[,"CTCF"]<-rbinom(N_PEAKS,1,0.5)
binding_temp[,"CTCF"]<-binding_status[,"CTCF"]
colnames(binding_status)<-FACTORS
colnames(binding_temp)<-FACTORS

# ---- Initialise parameters! ---- # These are all per-factor! Hence lists!
# this is also per-peak -> hence nested lists (ugh)
transition_params <- vector("list")
emission_params <- vector("list")
motifs<-vector("list")
for (factor in FACTORS){
    transition_params[[factor]] <- vector("list")
    motifs[[factor]]<-get_motif(factor)
}

#transition_params <- vector("list")
#for (factor in FACTORS){
#    transition_params[[factor]] <- vector("list")
#    for (peak in 1:N_PEAKS){
#        transition_params[[factor]][[peak]] <- list(B=c(0.4,0.2,0.4),G=c(0.5,0.5))
#    }
#}


# for testing: only looping over one TF
TEST_FACTORS<-"CTCF"
dbinding<-vector("numeric")
# ---- The outer loop: 'sample' over binding states ---- #
for (iter in 1:N_ITER){
    cat('Iteration',iter,'\n')
    # ---- Middle loop: iterate over each factor! ---- #
    for (factor in TEST_FACTORS){
        cat('Getting binding status for',factor,'\n')
        # get the motif
        pwm <- motifs[[factor]]
        factor_size <- ncol(pwm)
        N_STATES <- factor_size+2
        # initialise emission parameters for this factor
        if (!exists(paste("emission_params[[",factor,"]]",sep=""))){
            cat("No emission parameters saved for ",factor," - making some up!\n")
#            emission_params[[factor]] <- matrix(runif(N_STATES*N_FEATURES),nrow=N_STATES,ncol=N_FEATURES)
            emission_params[[factor]] <- matrix(rep(0.5,N_STATES*N_FEATURES),nrow=N_STATES,ncol=N_FEATURES)
        }
        # initialise a blank hmm with the right size... and fixed variables
        factor_hmm <- build_hmm(N_STATES,N_FEATURES,pwm)
        # calculate coincidence of this TF with the rest
        coincidence <- get_interactions(factor,binding_status)
        # when we finish EM we will record binding predictions for all peaks
        all_peaks_bound <- rep(NA,N_PEAKS)
        # initialise the while loop
        delta_ll <- EM_THRESHOLD*2
        ll_old <- -Inf
        EM.iter <- 0
        ll.all <- vector("numeric")
        decrease <- 0
        # ---- Second middle loop: EM! ---- #
        while(abs(delta_ll)>EM_THRESHOLD){
            EM.iter <- EM.iter + 1
            cat("EM iteration",EM.iter,"\n")
            theta_denom <- rep(0,N_STATES)
            theta_numer <- matrix(rep(0,N_STATES*N_FEATURES),nrow=N_STATES,ncol=N_FEATURES)
            # ll over the peaks
            ll_cumulative <- 0
            # ---- Inner loop: iterate over peaks! ---- #
            for (peak in 1:N_PEAKS){
                if (peak%%10000==0){
                    print(peak)
                }
                peak_data <- data[[peak]]
                peak_length <- ncol(peak_data)
               
                # initialise the parameters
                if (!exists(paste("transition_params[[",factor,"]][[",peak,"]]",sep=""))){
                    cat("No transmission parameters saved for ",factor," - making some up!\n")
                    transition_params[[factor]][[peak]] <- list(B=c(0.4,0.2,0.4),G=c(0.5,0.5))
                    }
                trans_param <- transition_params[[factor]][[peak]]
                emiss_param <- emission_params[[factor]]
                # set the parameters of the model
                peak_hmm <- initialise_hmm(factor_hmm,N_STATES,N_FEATURES,trans_param,emiss_param)
                # getting alpha and betas in here, basically
                gamma_and_xi <- get_gamma_and_xi(peak_hmm,peak_data,peak_length,N_STATES,N_FEATURES)
                theta <- get_theta(gamma_and_xi$"gamma",peak_data,N_STATES,N_FEATURES)
                #theta is additive over peaks, remember
                theta_numer <- theta_numer + theta$"theta_numer"
                theta_denom <- theta_denom + theta$"theta_denom"
#
                # get new transmission parameters (including a_BF)
                transition_params[[factor]][[peak]] <- get_new_transition_params(peak,peak_length,factor,factor_size,binding_status,coincidence,gamma_and_xi,TAU)

                # increase the log-likelihood...
                ll_cumulative <- ll_cumulative + gamma_and_xi$"ll"
            }
            # update emission parameters after all peaks
            # excluding emission_params for now... will EM converge?
#            emission_params[[factor]] <- theta_numer/theta_denom
#            emission_params[[factor]][is.nan(emission_params[[factor]])] <- 0
#            if(sum(emission_params[[factor]]==1)>0){
#                print('grr')
#                browser()
#            }

            # check how the likelihood has changed...
            ll <- ll_cumulative
            delta_ll <- ll-ll_old
            print(ll)
            if(delta_ll<0){
                cat("ERROR: likelihood is decreasing! Check yo EM!\n")
                print(delta_ll)
                decrease = decrease + delta_ll
            }
            ll_old <- ll
            ll.all <- c(ll.all,ll)
        }
        cat("EM has converged?\n")
        plot(ll.all,type='l',xlab='Iteration',ylab='Log-likelihood')
        cat('lhood decreased by',decrease,'in total\n')

        # now we have to check if it's bound or not...
        for (peak in 1:N_PEAKS){
            trans_param <- transition_params[[factor]][[peak]]
            emiss_param <- emission_params[[factor]]
            peak_hmm <- initialise_hmm(factor_hmm,N_STATES,N_FEATURES,trans_param,emiss_param)
            posteriors <- posterior.qhmm(peak_hmm,peak_data,n_threads=2)
            bound_yn <- is_bound(posteriors[,2])
            all_peaks_bound[peak] <- bound_yn
            }

        #path <- viterbi.qhmm(peak_hmm, peak_data)
        binding_temp[,factor] <- all_peaks_bound
    }
    # using the infinity norm here, for the... fun?
    dbinding <- c(dbinding,norm((binding_status-binding_temp),"I"))

    binding_status <- binding_temp

    # for the purpose of somehow gauging if convergence is occurring
#    visualise_binding(binding_status)
}
plot(dbinding)
# ---- After iteration: retrieve predictions ---- #
cm <- get_confusion_matrix(binding_status,bound_from_chip,FACTORS)
# This depends on how I'm storing the data, but basically need a prediction from each TF for each location, maybe whatever else... atm just doing it per-peak... can we do better than that? do we have a validation set?
