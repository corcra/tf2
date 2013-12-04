# Transcription factor binding site identification! Including interactions!
# NOTE THIS IS IN WEIRD PSEUDOCODE

# ---- Load functions! ---- #
source('wtf_fns.r')
library(rqhmm)

# ---- Load Data! ---- #
# DNase-I signal                            # does it make sense to have a peak class?
# DNA sequence
# PWM for TF(s)
#
# Ideas about peak class for DNAse-I signal!
# Do we need to do it peak by peak?
# Maybe maintain list of peak boundaries, then access the data as necessary.
# Array: LOC BASE [Dnase vars] ...?

# ---- Constants! ---- #
FACTORS <- c('f_one','f_two','f_three','f_four','f_five')
N_FACTORS <- length(FACTORS)
N_PEAKS <- 10
N_ITER <- 5
# how many features from DNAse will we take... (I'll be defining these somehow!)
N_FEATURES <- 4
all_the_data <- 100

# ---- Process the DNase-1 data! ---- #
#DNase_emissions <- get_features(DNase_data)

# ---- Initialise binding status ---- #
# Will need to get this from ChIP-seq data once I pick a test TF.
binding_status<-data.frame(matrix(rbinom(N_PEAKS*N_FACTORS,1,0.5),nrow=N_PEAKS,ncol=N_FACTORS))
colnames(binding_status)<-FACTORS

# ---- Initialise parameters! ---- #
initial_transition_params <- list(B=c(0.8,0.1,0.1),G=c(0.5,0.5))

# ---- The outer loop: 'sample' over binding states ---- #
for (iter in 1:N_ITER){
    cat('Iteration',iter,'\n')
    for (factor in FACTORS){
        cat('Getting binding status for',factor,'\n')
        pwm <- get_motif(factor)
        factor_size <- ncol(pwm)
        N_STATES <- factor_size+2
        # initialise a blank hmm with the right size... and fixed variables
        factor_hmm <- build_hmm(N_STATES,N_FEATURES,pwm)
        coincidence<-get_interactions(factor,binding_status)
        all_peaks_bound <- vector(length=N_PEAKS)
        for (peak in 1:N_PEAKS){
            cat('Peak',peak,'\n')
            peak_length <- 20
            a_BF <- get_a_BF(peak,peak_length,factor,binding_status,coincidence)
            cat('abf is:',a_BF,'\n')
            # build the HMM using Andre's library
            peak_hmm <- initialise_hmm(factor_hmm,a_BF,N_STATES,N_FEATURES,initial_transition_params)
            browser()
            # learn the other transitions with EM
            # FUNCTION BE HERE
            # this will be binary y/n
#            posteriors <- posterior.qhmm(peak_hmm,peak_data,n_threads=2)
            # on the basis of the posteriors, do we think the peak is bound?
            posteriors <- matrix(c(runif(peak_length*N_STATES)),nrow=peak_length,ncol=N_STATES)         # for now!
            all_peaks_bound[peak] <- is_bound(posteriors[,2])
            }
        #binding_status[[factor]]<-rbinom(N_PEAKS,1,0.5)
        binding_status[[factor]] <- all_peaks_bound
        }
    # for the purpose of somehow gauging if convergence is occurring
    visualise_binding(binding_status)
    }

# ---- After iteration: retrieve predictions ---- #
# This depends on how I'm storing the data, but basically need a prediction from each TF for each location, maybe whatever else.
