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

get_new_transition_params <- function(peak,peak_length,factor,factor_size,binding_status,coincidence,xis,TAU){
    # get the interaction modifier
    C_int <- get_C_int(peak,peak_length,factor,factor_size,binding_status,coincidence)

    norm_B <- xis$"xi_BB"+xis$"xi_BG"
    norm_G <- xis$"xi_GB"+xis$"xi_GG"

    a_BB <- (1-TAU-C_int)*xis$"xi_BB"/norm_B
    a_BF <- TAU + C_int
    a_BG <- (1-TAU-C_int)*xis$"xi_BG"/norm_B
    a_GB <- xis$"xi_GB"/norm_G
    a_GG <- xis$"xi_GG"/norm_G
   
    return(list(B=c(a_BB,a_BF,a_BG),G=c(a_GB,a_GG)))
}

get_theta <- function(gamma,data,N_STATES,N_FEATURES){
    theta_d <- rowSums(gamma)
    theta_n <- matrix(NA,nrow=N_STATES,ncol=N_FEATURES)
    # theta_numer is more complicated... for each location in the peak we have a matrix: the rows correspond to STATES and the columns correspond to the DNase emissions!
    for (feature in 1:N_FEATURES){
        # remember, ignoring the first emission (DNA sequence)
        X <- data[(feature+1),]-1
        for (state in 1:N_STATES){
            theta_n[state,feature]<-sum(gamma[state,]*X)
            if (is.na(theta_n[state,feature])){
                browser()
                }
            if (theta_n[state,feature]==theta_d[state]){
                # assume this is underflow... this is an ad hoc solution, i'm sorry :(
                theta_n[state,feature]<-0.95*theta_n[state,feature]+(0.05)*(1-theta_n[state,feature])
            }
        }
    }
# this is the 'smart' way of doing it... (without those middle 2 for loops)
#        theta_numer <- theta_numer + gamma[,loc]%o%(peak_data[2:(N_FEATURES+1),loc]-1)
    return(list("theta_numer"=theta_n,"theta_denom"=theta_d))
}

get_xi <- function(hmm,alpha_prime,beta_prime,loglik,peak_data,N_STATES,N_FEATURES){
    # get emission probabilities (this is the likelihood of the data given diff states, for each loc)
    emissions_B <- get_emission_prob(hmm,peak_data,peak_length,1,N_FEATURES)
    emissions_G <- get_emission_prob(hmm,peak_data,peak_length,N_STATES,N_FEATURES)
    # get transition probabilities
    a_BB <- get.transition.params.qhmm(hmm,1)[1]
    # remember, B can only transition to B,F,G
    a_BG <- get.transition.params.qhmm(hmm,1)[3]
    # G can only transition to B,G
    a_GB <- get.transition.params.qhmm(hmm,N_STATES)[1]
    a_GG <- get.transition.params.qhmm(hmm,N_STATES)[2]

    # the B state is the first row! G is final...
    xi_BB <- sum(exp(alpha_prime[1,1:(peak_length-1)]+beta_prime[1,2:peak_length]+log(emissions_B)+log(a_BB)-loglik))
    xi_BG <- sum(exp(alpha_prime[1,1:(peak_length-1)]+beta_prime[N_STATES,2:peak_length]+log(emissions_G)+log(a_BG)-loglik))
    xi_GB <- sum(exp(alpha_prime[N_STATES,1:(peak_length-1)]+beta_prime[1,2:peak_length]+log(emissions_B)+log(a_GB)-loglik))
    xi_GG <- sum(exp(alpha_prime[N_STATES,1:(peak_length-1)]+beta_prime[N_STATES,2:peak_length]+log(emissions_G)+log(a_GG)-loglik))

    return(list("xi_BB"=xi_BB,"xi_BG"=xi_BG,"xi_GB"=xi_GB,"xi_GG"=xi_GG))
}

get_confusion_matrix <- function(pred,known,FACTORS){
    cm <- vector("list")
    for (factor in FACTORS){
        pos<-which(known[,factor]==1)
        neg<-which(known[,factor]==0)
        TP<-sum(pred[pos,factor]==1)
        FP<-sum(pred[neg,factor]==1)
        TN<-sum(pred[neg,factor]==0)
        FN<-sum(pred[pos,factor]==0)
        cm[[factor]]<-list("TP"=TP,"FP"=FP,"TN"=TN,"FN"=FN,"sens"=100*TP/(TP+FN),"spec"=100*TN/(FP+TN))
        }
    return(cm)
}
