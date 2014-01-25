# Functions used in tf2.r!
cat('Loaded functions!\n')

# note! right now this only works for one factor...
get_motif<-function(factor){
    pfm <- read.table('data/CTCF.NGCGCCMYCTAGYGGTN.pwm',skip=2)
    # including some pseudo-counts... maybe unnecessary?
    pseudo_pfm <- matrix(rep(1,prod(dim(pfm))),dim(pfm))
    adj_pfm <- pfm+pseudo_pfm
    # note taking a transpose!
    pwm <- t(adj_pfm/rowSums(adj_pfm))
        return(pwm)
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

get_C_int<-function(peak,peak_length,factor,factor_size,binding_status,coincidence){
    # get the peak-specific bound status, excluding current factor
    peak_bound_status<-as.numeric(binding_status[peak,!colnames(binding_status)==factor])
    # translate this into a matrix of indices to query the coincidence matrix - the query will produce a vector of the relevant elements of the coincidence matrix depending on whether or not the factor is actully bound in this peak!
    relevant_indices<-matrix(c(seq(N_FACTORS-1),peak_bound_status+1),(N_FACTORS-1),2)
    # some normalisation, sum over aforementioned vector...
    C_int <- (1.0/((peak_length-factor_size)*(N_FACTORS-1)))*sum(coincidence[relevant_indices])
    return(C_int)
}

visualise_binding<-function(binding_status){
    print(binding_status)
}

build_hmm<-function(N_STATES,N_FEATURES,pwm){
    # no covariates!
    data_shape <- list(rep(1,N_FEATURES+1),NULL)

    # construct the valid transitions! this is a little tricky because the number of states depends on the transcription factor
    #valid_transitions <- matrix(rep(0,N_STATES*N_STATES),nrow=N_STATES,ncol=N_STATES,byrow=T)
    valid_transitions <- matrix(c(rep(1,N_STATES),rep(0,N_STATES*(N_STATES-1))),nrow=N_STATES,ncol=N_STATES,byrow=F)
    valid_transitions[1,] <- c(1,2,rep(0,(N_STATES-3)),3)
    #valid_transitions[2:(N_STATES-2),3:(N_STATES-1)] <- diag(N_STATES-3)
    valid_transitions[2:(N_STATES-2),3:(N_STATES-1)] <- 2*diag(N_STATES-3)
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
    #set.transition.params.qhmm(hmm,2:(N_STATES-1),1,fixed=T)
    set.transition.params.qhmm(hmm,2:(N_STATES-2),c(0.01,0.99),fixed=c(TRUE,TRUE))
    set.transition.params.qhmm(hmm,N_STATES-1,1,fixed=T)
    # include the pwm info!
    hmm <- include_pwm(hmm,pwm,N_STATES)
    return(hmm)
}

initialise_hmm<-function(hmm,N_STATES,N_FEATURES,transition_params,emission_params){
    # initialise it! ... always start in the background state!
    set.initial.probs.qhmm(hmm,c(1,rep(0,(N_STATES-1))))

    set.transition.params.qhmm(hmm,1,transition_params$B)
    set.transition.params.qhmm(hmm,N_STATES,transition_params$G)

    # emissions (DNAse)
    for (state in 1:N_STATES){
        for (feature in 1:N_FEATURES){
            theta <- emission_params[state,feature]
            if(theta==1|theta==0){
                print("we're about to get an emission probability of zero!")
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

# this just displays the likelihoods
visualise_probs <- function(gamma,peak,factor,N_STATES){
    rb<-rainbow(N_STATES-2)
    plot(gamma[1,],type="p",col="black",pch=1,cex=0.7,bty="l",ylab="Posterior",xlab="Location in peak")
    for(s in seq(2,N_STATES-1)){
        lines(gamma[s,],type="p",col=rb[s],pch=20)
    }
    lines(gamma[N_STATES,],type="p",col="grey",pch=20)
    browser()
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
                print(gamma[state,])
                print(X)
                browser()
                }
            if (theta_n[state,feature]==theta_d[state]){
                # assume this is underflow... this is an ad hoc solution, i'm sorry :(
                theta_n[state,feature]<-0.95*theta_n[state,feature]+(0.05)*(1-theta_n[state,feature])
            }
        }
    }
# this is the 'smart' way of doing it... (without those middle 2 for loops)
# theta_numer <- theta_numer + gamma[,loc]%o%(peak_data[2:(N_FEATURES+1),loc]-1)
    return(list("theta_numer"=theta_n,"theta_denom"=theta_d))
}

get_xi <- function(hmm,alpha_prime,beta_prime,loglik,peak_data,N_STATES,N_FEATURES,peak_length){
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

eval_peak <- function(peak,factor_hmm,data,trans_params,emiss_params,N_STATES,N_FEATURES,factor,factor_size,binding_status,coincidence,TAU){
    if (peak%%10000==0){
        print(peak)
    }

    peak_data <- data[[peak]]
    missing_data <- (peak_data==0)*1
    peak_length <- ncol(peak_data)

    # initialise parameters
    if (is.null(trans_params[[peak]])){
#        cat("No transition parameters saved for ",factor," - making some up!\n")
        trans_params[[peak]] <- list(B=c(0.4,0.2,0.4),G=c(0.5,0.5))
    }
    
    # initialise hmm
    if (sum(unlist(trans_params[[peak]])==0)>0){
        print("PROBLEM AHOY!")
        print(trans_params[[peak]])
    }

    peak_hmm <- initialise_hmm(factor_hmm,N_STATES,N_FEATURES,trans_params[[peak]],emiss_params)

    # get new parameters
    # rows are states, cols are locations
    alpha_prime <- forward.qhmm(peak_hmm,peak_data,missing=missing_data)
    beta_prime <- backward.qhmm(peak_hmm,peak_data,missing=missing_data)
    loglik <- attr(alpha_prime,"loglik")

    # emissions
    gamma <- exp(alpha_prime + beta_prime - loglik)
    theta <- get_theta(gamma,peak_data,N_STATES,N_FEATURES)

    # transitions
    xis <- get_xi(peak_hmm,alpha_prime,beta_prime,loglik,peak_data,N_STATES,N_FEATURES,peak_length)
    new_trans_params <- get_new_transition_params(peak,peak_length,factor,factor_size,binding_status,coincidence,xis,TAU)

    # check if bound (the result of this only makes sense after EM has converged)
    visualise_probs(gamma,peak,factor,N_STATES)
    posteriors <- gamma[2,]
    bound_yn <- is_bound(posteriors)
    
    # values to return
    return_vals <- list("theta"=theta,"new_trans_params"=new_trans_params,"loglik"=loglik,"bound"=bound_yn)
}
