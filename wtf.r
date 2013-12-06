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
N_PEAKS <- 206138
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
transition_params <- list(B=c(0.8,0.1,0.1),G=c(0.5,0.5))

# ---- The outer loop: 'sample' over binding states ---- #
for (iter in 1:N_ITER){
    cat('Iteration',iter,'\n')
# ---- Middle loop: iterate over each factor! ---- #
    for (factor in FACTORS){
        cat('Getting binding status for',factor,'\n')

        # get the motif
        pwm <- get_motif(factor)
        factor_size <- ncol(pwm)
        N_STATES <- factor_size+2

        # initialise a blank hmm with the right size... and fixed variables
        factor_hmm <- build_hmm(N_STATES,N_FEATURES,pwm)

        # calculate coincidence of this TF with the rest
        coincidence <- get_interactions(factor,binding_status)

        # when we finish EM we will record binding predictions for all peaks
        all_peaks_bound <- rep(NA,N_PEAKS)

# ---- Second middle loop: EM! ---- #
# THIS NEEDS SORTING OUT!
        ll <- 0
        theta_denom <- rep(0,N_STATES)
        # theta_numer is the numerator of theta, basically, gamma*X at each point... but there are N_FEATURES Xes remember! _-> hence this is a matrix _-> hence this is a matrix _-> hence this is a matrix _-> hence this is a matrix
        theta_numer <- matrix(rep(0,N_STATES*N_FEATURES),nrow=N_STATES,ncol=N_FEATURES)

        # testing running time...
#        print(system.time(for (peak in 1:N_PEAKS){get_a_BF(peak,peak_length,factor,factor_size,binding_status,coincidence)}))

# ---- Inner loop: iterate over peaks! ---- #
        for (peak in 1:N_PEAKS){
            # get the data
            # right now: fake data!
            peak_data <- matrix(as.numeric(rbinom((N_FEATURES+1)*100,1,0.5)),ncol=100,nrow=N_FEATURES+1)+1
            peak_length <- ncol(peak_data)

            # get a_BF!
            a_BF <- get_a_BF(peak,peak_length,factor,factor_size,binding_status,coincidence)
            #testing on known motif
#            motif <- c(1,1,1,3,2,3,2,2,1,2,2,4,1,3,4,3,3,4,1,1)
#            dnase_signal <- matrix(as.numeric(rbinom((N_FEATURES*20),1,0.5)),ncol=20,nrow=N_FEATURES)+1
#            peak_data <- rbind(motif,dnase_signal)
#            peak_length <- 20
#            cat('abf is:',a_BF,'\n')
#            if(is.na(a_BF)){
#                browser()}

            # build the HMM using Andre's library
            peak_hmm <- initialise_hmm(factor_hmm,a_BF,N_STATES,N_FEATURES,transition_params)
            
            # get the theta components (for DNase emissions, EM etc...)
            # also get the xis! (for transitions)
            # makes sense to do these at once because they both use the results from forward-backward
            theta_and_xi <- get_theta_and_xi(factor_hmm,peak_data,peak_length,N_STATES,N_FEATURES)

            theta_numer <- theta_numer + theta_and_xi$"theta_numer"
            theta_denom <- theta_denom + theta_and_xi$"theta_denom"

            # when I say xi, i really mean the summed form
            # translate these to transition probabilities!
            a_BB <- (1-a_BF)*theta_and_xi$"xi_BB"/(theta_and_xi$"xi_BB"+theta_and_xi$"xi_BG")
            a_BG <- (1-a_BF)*theta_and_xi$"xi_BG"/(theta_and_xi$"xi_BB"+theta_and_xi$"xi_BG")
            a_GB <- theta_and_xi$"xi_GB"/(theta_and_xi$"xi_GB"+theta_and_xi$"xi_GG")
            a_GG <- theta_and_xi$"xi_GG"/(theta_and_xi$"xi_GB"+theta_and_xi$"xi_GG")
            browser()

            # increase the log-likelihood...
            ll<-ll+theta_and_xi$"ll"

            posteriors <- posterior.qhmm(peak_hmm,peak_data,n_threads=2)
            #path <- viterbi.qhmm(peak_hmm, peak_data)
            browser()

            # on the basis of the posteriors, do we think the peak is bound? ... only do this after EM! ... iterate over peaks until convergence, then sweep through a final time to calculate the posteriors!
            #posteriors <- matrix(c(runif(peak_length*N_STATES)),nrow=peak_length,ncol=N_STATES)         # for now!
            bound_yn <- is_bound(posteriors[,2])
            if(is.na(bound_yn)){
                print("What's going on here?")
                print("It looks like I'm getting NaNs from the posterior call with this dataset... but why?")
                print("Think I figured it out... data encoding! It wasn't liking the zeroes.")
                browser()
                }
            all_peaks_bound[peak] <- bound_yn
            }
        #binding_status[[factor]]<-rbinom(N_PEAKS,1,0.5)
        binding_status[[factor]] <- all_peaks_bound
        }
    # for the purpose of somehow gauging if convergence is occurring
    }

visualise_binding(binding_status)
# ---- After iteration: retrieve predictions ---- #
# This depends on how I'm storing the data, but basically need a prediction from each TF for each location, maybe whatever else.
