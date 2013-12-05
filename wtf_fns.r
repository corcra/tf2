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
    hmm <- new.qhmm(data_shape,valid_transitions,transition_functions,emission_functions)

    # set the fixed parameters
    set.transition.params.qhmm(hmm,2:(N_STATES-1),1,fixed=T)
    # include the pwm info!
    hmm <- include_pwm(hmm,pwm,N_STATES)
    return(hmm)
    }

initialise_hmm<-function(hmm,a_BF,N_STATES,N_FEATURES,initial_transition_params){
    # initialise it!
    set.initial.probs.qhmm(hmm,c(1,rep(0,(N_STATES-1))))

    # transitions
    set.transition.params.qhmm(hmm,1,initial_transition_params$B)
    set.transition.params.qhmm(hmm,N_STATES,initial_transition_params$G)

    # emissions (DNAse)
    # NOTE! WHAT'S GOING ON HERE? do we have any idea of prior? can we learn this as well? should we learn this as well?
    set.emission.params.qhmm(hmm, 1:N_STATES, rep(0.5,2), slot=2:(N_FEATURES+1))
    return(hmm)
    }


# coincidence matrix is one row per transcription factor, first column is 'F's CDF given no phi', second is 'F's CDF with phi' approximately speaking...
# since we only update one factor at a time, it's wasteful to recalculate for all the others every time, but I'll worry about that some other time
get_interactions<-function(factor,binding_status){
    # binding of factor of interest
    reference<-binding_status[[factor]]
    # all the others
    rest_mat<-binding_status[,!names(binding_status)==factor]
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
    coincidence_mat[is.nan(coincidence_mat)]=0
    return(coincidence_mat)
    }

get_a_BF<-function(peak,peak_length,factor,binding_status,coincidence){
    # get the peak-specific bound status, excluding current factor
    peak_bound_status<-as.numeric(binding_status[peak,!names(binding_status)==factor])
    # translate this into a matrix of indices to query the coincidence matrix - the query will produce a vector of the relevant elements of the coincidence matrix depending on whether or not the factor is actully bound in this peak!
    relevant_indices<-matrix(c(seq(N_FACTORS-1),peak_bound_status+1),(N_FACTORS-1),2)
    # some normalisation, sum over aforementioned vector...
    a_BF <- (1.0/((peak_length-factor_size)*(N_FACTORS-1)))*sum(coincidence[relevant_indices])
    return(a_BF)
    }

visualise_binding<-function(binding_status){
    print(binding_status)
    }
