# Transcription factor binding site identification! Including interactions!

# ---- Constants! ---- #
FACTORS <- c('f_one','f_two','f_three','f_four','f_five')
N_FACTORS <- length(FACTORS)
#N_PEAKS <- 112025
N_PEAKS <- 20
N_ITER <- 5
# DNase features
N_FEATURES <- 1
THRESHOLD <- 0.5
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
# Question: how much preprocessing to do on the data? Definitely need:
#   - sequence (in FASTA format)
#   - DNase values! (eg a bw or wiggle file!)
# Or maybe just take DNase peak calls (+wig file), get the seq from that...
#
# right now this data is real sequence data but made-up features
fc <- file('processed_data.gz',open='r')
data <- vector("list",N_PEAKS)
for (peak in 1:N_PEAKS){
    buff <- scan(fc,sep=" ",what=numeric(),nlines=(N_FEATURES+1))
    # columns -> number of locations, rows -> number of emission variables (first one will be DNA)
    data[[peak]] <- matrix(buff,nrow=(N_FEATURES+1),byrow=T)
    len<-ncol(data[[peak]])
    #data[[peak]][2,] <- rbinom(len,1,0.5)+1
    data[[peak]][2,] <- seq(len)%%2+1
    }
close(fc)
cat("Data loaded!\n")

# ---- Initialise binding status ---- #
# Will need to get this from ChIP-seq data once I pick a test TF.
binding_status<-matrix(rbinom(N_PEAKS*N_FACTORS,1,0.5),nrow=N_PEAKS,ncol=N_FACTORS)
binding_temp<-matrix(rbinom(N_PEAKS*N_FACTORS,1,0.5),nrow=N_PEAKS,ncol=N_FACTORS)
#binding_temp<-matrix(rep(0,(N_PEAKS*N_FACTORS)),nrow=N_PEAKS,ncol=N_FACTORS)
colnames(binding_status)<-FACTORS
colnames(binding_temp)<-FACTORS

# ---- Initialise parameters! ---- # These are all per-factor! Hence lists!
# this is also per-peak -> hence nested lists (ugh)
transition_params <- vector("list")
for (factor in FACTORS){
    transition_params[[factor]] <- vector("list")
    for (peak in 1:N_PEAKS){
        transition_params[[factor]][[peak]] <- list(B=c(0.4,0.2,0.4),G=c(0.5,0.5))
    }
}

emission_params <- vector("list")

motifs<-vector("list")
for (factor in FACTORS){
    motifs[[factor]]<-get_motif(factor)
}

# for testing: only looping over one TF
TEST_FACTORS<-"f_one"
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
            print(paste("no emission parameters saved for",factor,"- making some up!"))
            emission_params[[factor]] <- matrix(runif(N_STATES*N_FEATURES),nrow=N_STATES,ncol=N_FEATURES)
            }

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
        ll.all <- vector("numeric")

        decrease <- 0
        # ---- Second middle loop: EM! ---- #
        while(abs(delta_ll)>THRESHOLD){
            EM.iter <- EM.iter + 1
            cat("EM iteration",EM.iter,"\n")
            theta_denom <- rep(0,N_STATES)
            # theta_numer is the numerator of theta, basically, gamma*X at each point... but there are N_FEATURES Xes remember! _-> hence this is a matrix _-> hence this is a matrix _-> hence this is a matrix _-> hence this is a matrix
            theta_numer <- matrix(rep(0,N_STATES*N_FEATURES),nrow=N_STATES,ncol=N_FEATURES)

            # ll over the peaks
            ll_cumulative <- 0
            test_gamma <- rep(0,N_STATES)
            # ---- Inner loop: iterate over peaks! ---- #
            for (peak in 1:N_PEAKS){
                if (peak%%10000==0){
                    print(peak)
                }
                peak_data <- data[[peak]]
                if(mean(peak_data[2,]==1)==1){
                    print("all 1s...")
                    browser()
                    }
                if(mean(peak_data[2,]==2)==1){
                    print("all 2s...")
                    browser()
                    }
 
                peak_length <- ncol(peak_data)
               
                # initialise the parameters
                trans_param <- transition_params[[factor]][[peak]]
                emiss_param <- emission_params[[factor]]
                peak_hmm <- initialise_hmm(factor_hmm,N_STATES,N_FEATURES,trans_param,emiss_param)
                
                # get the theta components (for DNase emissions, EM etc...)
                # also get the xis! (for transitions)
                # makes sense to do these at once because they both use the results from forward-backward
                gamma_and_xi <- get_gamma_and_xi(factor_hmm,peak_data,peak_length,N_STATES,N_FEATURES)

                test_gamma <- test_gamma + rowSums(gamma_and_xi$"gamma")
                theta <- get_theta(gamma_and_xi$"gamma",peak_data,N_STATES,N_FEATURES)
                # incease the theta counts ... will collect all of these at the end of the peak loop
                theta_numer <- theta_numer + theta$"theta_numer"
                theta_denom <- theta_denom + theta$"theta_denom"
    
#                if(sum(theta_numer==theta_denom)>0){
#                    print('theta_num is the same as theta_denom somewhere?')
#                    browser()
#                }
#
                # save the transition parameters for this peak (we will use these next time)
                # note! this includes getting the interaction part!
                transition_params[[factor]][[peak]] <- get_new_transition_params(peak,peak_length,factor,factor_size,binding_status,coincidence,gamma_and_xi,TAU)

                # increase the log-likelihood...
                ll_cumulative <- ll_cumulative + gamma_and_xi$"ll"
            }
            # update the emission parameters ... the transition parameters are saved in transition_params
#            if(sum(theta_denom==0)>0){
#                print("have zeroes in theta_denom!")
#                browser()
#            }
#            if(sum(theta_numer==0)>0){
#                print("have zeroes in theta_numer!")
#                browser()
#            }
# excluding emission_params for now... will EM converge? who knows!
#            for (state in 1:N_STATES){
#                emission_params[[factor]][state,] <- theta_numer[state]/theta_denom[state]
#            }
# 'smart' way of doing this
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
#                browser()
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
#            print(bound_yn)
            }

        #path <- viterbi.qhmm(peak_hmm, peak_data)
        binding_temp[,factor] <- all_peaks_bound
    }
    binding_status <- binding_temp
    # for the purpose of somehow gauging if convergence is occurring
    visualise_binding(binding_status)
}
# ---- After iteration: retrieve predictions ---- #
# This depends on how I'm storing the data, but basically need a prediction from each TF for each location, maybe whatever else.
