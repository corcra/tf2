# This code contains functions for wtf.r
cat('Loaded potentially-dodgy functions!\n')

get_features<-function(DNase_data){
# what are the features, eh?
}

build_hmm<-function(N_STATES,N_FEATURES,pwm){
# no covariates!
    data_shape <- list(rep(1,N_FEATURES+1),NULL)

# construct the valid transitions! this is a little tricky because the number of states depends on the transcription factor
    valid_transitions <- matrix(rep(0,N_STATES*N_STATES),nrow=N_STATES,ncol=N_STATES,byrow=T)
    valid_transitions[1,] <- c(1,2,rep(0,(N_STATES-3)),3)
    valid_transitions[2:(N_STATES-2),3:(N_STATES-1)] <- diag(N_STATES-3)
    valid_transitions[(N_STATES-1),] <- c(1,rep(0,(N_STATES-1)))
    valid_transitions[N_STATES,] <- c(1, rep(0,(N_STATES-2)), 2)
# everything is discrete here!
    transition_functions <- rep("discrete",N_STATES)
# noo, a for loop!
    emission_functions <- vector(mode="list",length=N_STATES)
    for (state in 1:N_STATES){
        emission_functions[[state]]<-rep("discrete",(N_FEATURES+1))
    }
# create the HMM
    hmm <- new.qhmm(data_shape,valid_transitions,transition_functions,emission_functions,support.missing=TRUE)

# set the fixed parameters
    set.transition.params.qhmm(hmm,2:(N_STATES-1),1,fixed=T)
# include the pwm info!
    hmm <- include_pwm(hmm,pwm,N_STATES)
    return(hmm)
}

initialise_hmm<-function(hmm,N_STATES,N_FEATURES,transition_params,emission_params){
    # initialise it! ... always start in the background state!
    set.initial.probs.qhmm(hmm,c(1,rep(0,(N_STATES-1))))

    # transitions ... the only non-fixed ones are from B and G, we give these as input!
    if (sum(is.nan(transition_params$G))==0){
        set.transition.params.qhmm(hmm,1,transition_params$B)
        set.transition.params.qhmm(hmm,N_STATES,transition_params$G)
    }
    else{
        browser()
    }

    # emissions (DNAse)
    for (state in 1:N_STATES){
        for (feature in 1:N_FEATURES){
            theta <- emission_params[state,feature]
            if(theta==0|theta==1){
                print("we're about to get an emission probability of zero...")
                }
            set.emission.params.qhmm(hmm, state, c(1-theta,theta),slot=(feature+1))
            }
        }
    return(hmm)
}

get_emission_prob<-function(hmm,peak_data,peak_length,state,N_FEATURES){
    probs<-rep(1,(peak_length-1))
    for (slot_num in 1:(N_FEATURES+1)){
        emission_probs<-get.emission.params.qhmm(hmm,state,slot=slot_num)
        # this line might look a bit weird, but the value in the peak_data is actually the label found in the hmm
        # also note: we're actually only interested in the probs from 2... end... see the mathematical form!
        probs<-probs*emission_probs[peak_data[slot_num,2:peak_length]]
    }
    return(probs)
}

get_new_transition_params <- function(peak,peak_length,factor,factor_size,binding_status,coincidence,theta_and_xi){
    # get the interaction modifier
    C_int <- get_C_int(peak,peak_length,factor,factor_size,binding_status,coincidence)

    # when I say xi, i really mean the summed form
    # translate these to transition probabilities!
    norm_B <- theta_and_xi$"xi_BB"+theta_and_xi$"xi_BF"+theta_and_xi$"xi_BG"
    norm_G <- theta_and_xi$"xi_GB"+theta_and_xi$"xi_GG"
    if (norm_B==0|norm_G==0){
        print('wtf')
        browser()
        }
    a_BB <- (1-C_int)*theta_and_xi$"xi_BB"/norm_B
    a_BF <- (1-C_int)*theta_and_xi$"xi_BF"/norm_B + C_int
    a_BG <- (1-C_int)*theta_and_xi$"xi_BG"/norm_B
    a_GB <- theta_and_xi$"xi_GB"/norm_G
    a_GG <- theta_and_xi$"xi_GG"/norm_G
   
    if(!sum(is.na(c(a_BB,a_BG,a_GB,a_GG)))==0){
        browser()
        }
    return(list(B=c(a_BB,a_BF,a_BG),G=c(a_GB,a_GG)))
}

get_theta_and_xi<-function(hmm,peak_data,peak_length,N_STATES,N_FEATURES){
    # each ROW of this corresponds to a different state
    alpha_prime <- forward.qhmm(factor_hmm,peak_data)
    beta_prime <- backward.qhmm(factor_hmm,peak_data)
    ll <- attr(alpha_prime,"loglik")
    # if it's not finite it's probably -Inf, which is secretly a 0
    alpha_prime[!is.finite(alpha_prime)]<-0
    beta_prime[!is.finite(beta_prime)]<-0

    # THETA !
    # this here is a vector! ... remember the running total
    gamma <- exp(alpha_prime + beta_prime-ll)
    theta_denom <- rowSums(gamma)

    # theta_numer is more complicated... for each location in the peak we have a matrix: the rows correspond to STATES and the columns correspond to the DNase emissions!
    # another cursed for loop!
    theta_numer <- 0
    for (loc in 1:peak_length){
    # note: the need the emissions to be truly 0 and 1 here, not 1 and 2 as we gave it to the qhmm!
        theta_numer <- theta_numer + gamma[,loc]%o%(peak_data[2:(N_FEATURES+1),loc]-1)
        }

    # XI !
    # Note: only care about xi_BB, (xi_BF), xi_BG, xi_GB, and xi_GG, so can do these separately...
    # Will be needing the likelihood of the data given G and B...
    emissions_B <- get_emission_prob(hmm,peak_data,peak_length,1,N_FEATURES)
    emissions_F <- get_emission_prob(hmm,peak_data,peak_length,2,N_FEATURES)
    emissions_G <- get_emission_prob(hmm,peak_data,peak_length,N_STATES,N_FEATURES)
    a_BB <- get.transition.params.qhmm(hmm,1)[1]
    a_BF <- get.transition.params.qhmm(hmm,1)[2]
    # remember, B can only transition to B,F,G
    a_BG <- get.transition.params.qhmm(hmm,1)[3]
    # G can only transition to B,G
    a_GB <- get.transition.params.qhmm(hmm,N_STATES)[1]
    a_GG <- get.transition.params.qhmm(hmm,N_STATES)[2]

    # the B state is the first row!
    xi_BB <- sum(exp(alpha_prime[1,1:(peak_length-1)]+beta_prime[1,2:peak_length]+log(emissions_B)+log(a_BB)-ll))
    xi_BF <- sum(exp(alpha_prime[1,1:(peak_length-1)]+beta_prime[2,2:peak_length]+log(emissions_F)+log(a_BF)-ll))
    xi_BG <- sum(exp(alpha_prime[1,1:(peak_length-1)]+beta_prime[N_STATES,2:peak_length]+log(emissions_G)+log(a_BG)-ll))
    xi_GB <- sum(exp(alpha_prime[N_STATES,1:(peak_length-1)]+beta_prime[1,2:peak_length]+log(emissions_B)+log(a_GB)-ll))
    xi_GG <- sum(exp(alpha_prime[N_STATES,1:(peak_length-1)]+beta_prime[N_STATES,2:peak_length]+log(emissions_G)+log(a_GG)-ll))

    return(list("theta_numer"=theta_numer,"theta_denom"=theta_denom,"xi_BB"=xi_BB,"xi_BG"=xi_BG,"xi_BF"=xi_BF,"xi_GB"=xi_GB,"xi_GG"=xi_GG,"ll"=ll))
}
