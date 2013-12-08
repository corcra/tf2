# This code contains functions-under-development for wtf.r
cat('Loaded potentially-dodgy functions!\n')

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

get_new_transition_params <- function(peak,peak_length,factor,factor_size,binding_status,coincidence,gamma_and_xi){
    # get the interaction modifier
    C_int <- get_C_int(peak,peak_length,factor,factor_size,binding_status,coincidence)

    # when I say xi, i really mean the summed form
    # translate these to transition probabilities!
    norm_B <- gamma_and_xi$"xi_BB"+gamma_and_xi$"xi_BF"+gamma_and_xi$"xi_BG"
    norm_G <- gamma_and_xi$"xi_GB"+gamma_and_xi$"xi_GG"
    if (norm_B==0|norm_G==0){
        print('wtf')
        browser()
        }
    a_BB <- (1-C_int)*gamma_and_xi$"xi_BB"/norm_B
    a_BF <- (1-C_int)*gamma_and_xi$"xi_BF"/norm_B + C_int
    a_BG <- (1-C_int)*gamma_and_xi$"xi_BG"/norm_B
    a_GB <- gamma_and_xi$"xi_GB"/norm_G
    a_GG <- gamma_and_xi$"xi_GG"/norm_G
   
    if(!sum(is.na(c(a_BB,a_BG,a_GB,a_GG)))==0){
        browser()
        }
    return(list(B=c(a_BB,a_BF,a_BG),G=c(a_GB,a_GG)))
}

get_theta <- function(gamma,peak_data,N_STATES,N_FEATURES){
    theta_d <- rowSums(gamma)
    theta_n <- matrix(rep(0,N_STATES*N_FEATURES),nrow=N_STATES,ncol=N_FEATURES)
    # theta_numer is more complicated... for each location in the peak we have a matrix: the rows correspond to STATES and the columns correspond to the DNase emissions!
    
    # another cursed for loop!
    for (loc in 1:peak_length){
    # note: the need the emissions to be truly 0 and 1 here, not 1 and 2 as we gave it to the qhmm!
        for (state in 1:N_STATES){
            for (feature in 1:N_FEATURES){
#                print('vals:')
#                print(gamma[state,loc])
#                print(peak_data[(feature+1),loc]-1)
#                print('before:')
#                print(theta_n[state,feature])
                theta_n[state,feature] <- theta_n[state,feature] + gamma[state,loc]*(peak_data[(feature+1),loc]-1)
#                print("after:")
#                print(theta_n[state,feature])
#                print("because:")
#                print(gamma[state,loc]*(peak_data[(feature+1),loc]-1))
#                browser()
            }
        }
    }
# this is the 'smart' way of doing it... (without those middle 2 for loops)
#        theta_numer <- theta_numer + gamma[,loc]%o%(peak_data[2:(N_FEATURES+1),loc]-1)
    return(list("theta_numer"=theta_n,"theta_denom"=theta_d))
}

get_gamma_and_xi<-function(hmm,peak_data,peak_length,N_STATES,N_FEATURES){
    # each ROW of this corresponds to a different state
    alpha_prime <- forward.qhmm(factor_hmm,peak_data)
    beta_prime <- backward.qhmm(factor_hmm,peak_data)
    loglik <- attr(alpha_prime,"loglik")
#    # if it's not finite it's probably -Inf, which is secretly a 0
#    wait, we don't want ot delete these...
#    alpha_prime[!is.finite(alpha_prime)]<-0
#    beta_prime[!is.finite(beta_prime)]<-0

    # THETA !
    # this here is a vector! ... remember the running total
    gamma <- exp(alpha_prime + beta_prime-loglik)
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
    xi_BB <- sum(exp(alpha_prime[1,1:(peak_length-1)]+beta_prime[1,2:peak_length]+log(emissions_B)+log(a_BB)-loglik))
    xi_BF <- sum(exp(alpha_prime[1,1:(peak_length-1)]+beta_prime[2,2:peak_length]+log(emissions_F)+log(a_BF)-loglik))
    xi_BG <- sum(exp(alpha_prime[1,1:(peak_length-1)]+beta_prime[N_STATES,2:peak_length]+log(emissions_G)+log(a_BG)-loglik))
    xi_GB <- sum(exp(alpha_prime[N_STATES,1:(peak_length-1)]+beta_prime[1,2:peak_length]+log(emissions_B)+log(a_GB)-loglik))
    xi_GG <- sum(exp(alpha_prime[N_STATES,1:(peak_length-1)]+beta_prime[N_STATES,2:peak_length]+log(emissions_G)+log(a_GG)-loglik))

    return(list("gamma"=gamma,"xi_BB"=xi_BB,"xi_BG"=xi_BG,"xi_BF"=xi_BF,"xi_GB"=xi_GB,"xi_GG"=xi_GG,"ll"=loglik))
}
