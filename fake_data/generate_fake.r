# Some functions to make it easier to generate fake data!

bases<-function(n){ return(ceiling(4*runif(n))) }

signal<-function(n){
    s1<-rbinom(n,1,0.5)+1
    s<-rbind(s1,rbinom(n,1,0.5)+1)
    s<-rbind(s,rbinom(n,1,0.5)+1)
    return(s)
}

save_df<-function(df,name){ write.table(df,file=name,col.names=F,row.names=F,append=T)
}

mutate<-function(motif,r){
    n<-length(motif)
    rs<-runif(n)
    new<-bases(n)
    change<-which(rs<r)
    motif[change]<-new[change]
    return(motif)
}

# note this is the CTCF motif
# 1,2,3,4: A,C,G,T
motif<-c(1, 3, 2, 3, 2, 2, 1, 2, 2, 4, 1, 3, 4, 3, 3, 4, 1)

with_mutations<-c(bases(5),mutate(motif,0),bases(5),mutate(motif,0.1),bases(5),mutate(motif,0.2),bases(5),mutate(motif,0.3),bases(5),mutate(motif,0.4),bases(5),mutate(motif,0.5),bases(5),mutate(motif,0.6),bases(5),mutate(motif,0.7),bases(5),mutate(motif,0.8),bases(5),mutate(motif,0.9),bases(5))
df_mutated<-rbind(with_mutations,signal(length(with_mutations)))

