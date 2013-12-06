# Transcription factor binding site identification! Including interactions!

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
N_PEAKS <- 20
N_ITER <- 5
# how many features from DNAse will we take... (I'll be defining these somehow!)
N_FEATURES <- 4
THRESHOLD <- 0.01

# ---- Process the DNase-1 data! ---- #
#DNase_emissions <- get_features(DNase_data)

# ---- Initialise binding status ---- #
# Will need to get this from ChIP-seq data once I pick a test TF.
binding_status<-matrix(rbinom(N_PEAKS*N_FACTORS,1,0.5),nrow=N_PEAKS,ncol=N_FACTORS)
binding_temp<-matrix(rep(0,(N_PEAKS*N_FACTORS)),nrow=N_PEAKS,ncol=N_FACTORS)
colnames(binding_status)<-FACTORS
colnames(binding_temp)<-FACTORS

# ---- Initialise parameters! ---- #
# note here: this is a list... one set of parameters for each peak!
transition_params <- vector("list",N_PEAKS)
# set up initial values
for (peak in 1:N_PEAKS){
    transition_params[[peak]]<- list(B=c(0.8,0.1,0.1),G=c(0.5,0.5))
    }
# set up initial emission params... what should I choose here?
emission_params <- matrix(runif(N_STATES*N_FEATURES),nrow=N_STATES,ncol=N_FEATURES)
# record the a_BFs for each peak, too
all_a_BFs <- vector("numeric",N_PEAKS)

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

        # initialise the while loop
        delta_ll <- THRESHOLD*2
        ll_old <- -Inf
        EM.iter <- 0

# ---- Second middle loop: EM! ---- #
        while(delta_ll>THRESHOLD){
            EM.iter <- EM.iter + 1
            cat("Iter",EM.iter,"\n")
            theta_denom <- rep(0,N_STATES)
            # theta_numer is the numerator of theta, basically, gamma*X at each point... but there are N_FEATURES Xes remember! _-> hence this is a matrix _-> hence this is a matrix _-> hence this is a matrix _-> hence this is a matrix
            theta_numer <- matrix(rep(0,N_STATES*N_FEATURES),nrow=N_STATES,ncol=N_FEATURES)

            # ll over the peaks
            ll_cumulative <- 0
    # ---- Inner loop: iterate over peaks! ---- #
            for (peak in 1:N_PEAKS){
                # get the data
                peak_data <- matrix(as.numeric(rbinom((N_FEATURES+1)*100,1,0.5)),ncol=100,nrow=N_FEATURES+1)+1
                peak_length <- ncol(peak_data)
               
                # initialise the parameters
                #peak_hmm <- initialise_hmm(factor_hmm,N_STATES,N_FEATURES,transition_params[[peak]],emission_params)
                # obvs fake
                peak_hmm <- initialise_hmm(factor_hmm,N_STATES,N_FEATURES,transition_params[[peak]],emission_params)
                
                # get the theta components (for DNase emissions, EM etc...)
                # also get the xis! (for transitions)
                # makes sense to do these at once because they both use the results from forward-backward
                theta_and_xi <- get_theta_and_xi(factor_hmm,peak_data,peak_length,N_STATES,N_FEATURES)

                # get a_BF!
                a_BF <- get_a_BF(peak,peak_length,factor,factor_size,binding_status,coincidence)
 
                # save the transition parameters for this peak (we will use these next time)
                transition_params[[peak]] <- get_new_transition_params(theta_and_xi,a_BF)

                # incease the theta counts ... will collect all of these at the end of the peak loop
                theta_numer <- theta_numer + theta_and_xi$"theta_numer"
                theta_denom <- theta_denom + theta_and_xi$"theta_denom"

                # increase the log-likelihood...
                ll_cumulative <- ll_cumulative +theta_and_xi$"ll"
                }
            # update the emission parameters ... the transition parameters are saved in transition_params
            emission_params <- theta_numer/theta_denom

            # check how the likelihood has changed...
            ll <- ll_cumulative
            delta_ll <- ll-ll_old
            print(ll)
            if(delta_ll<0){
                cat("ERROR: likelihood is decreasing! Check yo EM!\n")
                browser()
                }
            ll_old <- ll
            }
            cat("EM has converged?\n")
            # by the time we get here, EM has converged ... hopefully!

 #            posteriors <- posterior.qhmm(peak_hmm,peak_data,n_threads=2)
            #path <- viterbi.qhmm(peak_hmm, peak_data)
 #           browser()

            # on the basis of the posteriors, do we think the peak is bound? ... only do this after EM! ... iterate over peaks until convergence, then sweep through a final time to calculate the posteriors!
            #posteriors <- matrix(c(runif(peak_length*N_STATES)),nrow=peak_length,ncol=N_STATES)         # for now!
 #           bound_yn <- is_bound(posteriors[,2])
 #           all_peaks_bound[peak] <- bound_yn
            
        binding_temp[,factor] <- all_peaks_bound
        }
    binding_status <- binding_temp
    # for the purpose of somehow gauging if convergence is occurring
    visualise_binding(binding_status)
    }

# ---- After iteration: retrieve predictions ---- #
# This depends on how I'm storing the data, but basically need a prediction from each TF for each location, maybe whatever else.
