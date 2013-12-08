# This code contains functions for wtf.r
cat('Loaded functions!\n')

# note! right now this only works for one factor... but could be extended
get_motif<-function(factor){
    pfm <- read.table('CTCF.NGCGCCMYCTAGYGGTN.pwm',skip=2)
# including some pseudo-counts... maybe unnecessary?
        pseudo_pfm <- matrix(rep(1,prod(dim(pfm))),dim(pfm))
        adj_pfm <- pfm+pseudo_pfm
# note taking a transpose!
        pwm <- t(adj_pfm/rowSums(adj_pfm))
        return(pwm)
}

get_features<-function(DNase_data){
# what are the features, eh?
}

is_bound<-function(probs){
    success <- runif(length(probs))<probs
# a single success is sufficient for binding!
        bound_yn <- ifelse((sum(success)==0),0,1)
        return(bound_yn)
}

include_pwm<-function(hmm,pwm,N_STATES){
# assume background and generic states have equal probs
    set.emission.params.qhmm(hmm, c(1,N_STATES), rep(0.25,4), slot=1, fixed=rep(T,4))
# the pwm is just a matrix, whose columns refer to positions and rows refer to ACGT
        for (state in 2:(N_STATES-1)){
            set.emission.params.qhmm(hmm, state, pwm[,(state-1)], slot=1, fixed=rep(T,4))
        }
    return(hmm)
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
    set.transition.params.qhmm(hmm,1,transition_params$B)
    set.transition.params.qhmm(hmm,N_STATES,transition_params$G)

    # emissions (DNAse)
    for (state in 1:N_STATES){
        for (feature in 1:N_FEATURES){
            theta <- emission_params[state,feature]
            set.emission.params.qhmm(hmm, state, c(1-theta,theta),slot=(feature+1))
            }
        }
    return(hmm)
}


# coincidence matrix is one row per transcription factor, first column is 'F's CDF given no phi', second is 'F's CDF with phi' approximately speaking...
get_interactions<-function(factor,binding_status){
# binding of factor of interest
    reference<-binding_status[,factor]
# all the others
        rest_mat<-binding_status[,!colnames(binding_status)==factor]
# if they're both bound, we see 2!
        comp_mat<-rest_mat+reference
# if other isn't bound, but ref is, we see -1!
        cont_mat<-rest_mat-reference
# how often do we see both, given we saw 'other'?
        if_yes<-colSums(comp_mat==2)/colSums(rest_mat==1)
# how often do we see ref and not other, given no other?
        if_no<-colSums(cont_mat==-1)/colSums(rest_mat==0)
# turn this into a matrix, and pass it back
        coincidence_mat<-matrix(c(if_no,if_yes),nrow=(N_FACTORS-1),ncol=2)
# remove the NaNs... they should be zeroes anyway!
        coincidence_mat[is.nan(coincidence_mat)]<-0
        return(coincidence_mat)
}

get_a_BF<-function(peak,peak_length,factor,factor_size,binding_status,coincidence){
# get the peak-specific bound status, excluding current factor
    peak_bound_status<-as.numeric(binding_status[peak,!colnames(binding_status)==factor])
# translate this into a matrix of indices to query the coincidence matrix - the query will produce a vector of the relevant elements of the coincidence matrix depending on whether or not the factor is actully bound in this peak!
        relevant_indices<-matrix(c(seq(N_FACTORS-1),peak_bound_status+1),(N_FACTORS-1),2)
# some normalisation, sum over aforementioned vector...
        a_BF <- (1.0/((peak_length-factor_size)*(N_FACTORS-1)))*sum(coincidence[relevant_indices])
        return(a_BF)
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

get_new_transition_params <- function(theta_and_xi,a_BF){
    # when I say xi, i really mean the summed form
    # translate these to transition probabilities!
    norm_B <- theta_and_xi$"xi_BB"+theta_and_xi$"xi"
    norm_G <- theta_and_xi$"xi_GB"+theta_and_xi$"xi_GG"
    if (norm_B>0){
        a_BB <- (1-a_BF)*theta_and_xi$"xi_BB"/norm_B
        a_BG <- (1-a_BF)*theta_and_xi$"xi_BG"/norm_B
    }
    else{
        a_BB <- 0
        a_BG <- 0
    }
    if (norm_G>0){
        a_GB <- theta_and_xi$"xi_GB"/norm_G
        a_GG <- theta_and_xi$"xi_GG"/norm_G
    }
    else{
        a_GB <- 0
        a_GG <- 0
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
    # Note: only care about xi_BB, xi_BG, xi_GB, and xi_GG, so can do these separately...
    # Will be needing the likelihood of the data given G and B...
    emissions_B <- get_emission_prob(hmm,peak_data,peak_length,1,N_FEATURES)
    emissions_G <- get_emission_prob(hmm,peak_data,peak_length,N_STATES,N_FEATURES)
    a_BB <- get.transition.params.qhmm(hmm,1)[1]
    # remember, B can only transition to B,F,G
    a_BG <- get.transition.params.qhmm(hmm,1)[3]
    # G can only transition to B,G
    a_GG <- get.transition.params.qhmm(hmm,N_STATES)[2]
    a_GB <- get.transition.params.qhmm(hmm,N_STATES)[1]

    # the B state is the first row!
    xi_BB <- sum(exp(alpha_prime[1,1:(peak_length-1)]+beta_prime[1,2:peak_length]+log(emissions_B)+log(a_BB)-ll))
    xi_BG <- sum(exp(alpha_prime[1,1:(peak_length-1)]+beta_prime[N_STATES,2:peak_length]+log(emissions_G)+log(a_BG)-ll))
    xi_GB <- sum(exp(alpha_prime[N_STATES,1:(peak_length-1)]+beta_prime[1,2:peak_length]+log(emissions_B)+log(a_GB)-ll))
    xi_GG <- sum(exp(alpha_prime[N_STATES,1:(peak_length-1)]+beta_prime[N_STATES,2:peak_length]+log(emissions_G)+log(a_GG)-ll))

    return(list("theta_numer"=theta_numer,"theta_denom"=theta_denom,"xi_BB"=xi_BB,"xi_BG"=xi_BG,"xi_GB"=xi_GB,"xi_GG"=xi_GG,"ll"=ll))
}

visualise_binding<-function(binding_status){
    print(binding_status)
}
