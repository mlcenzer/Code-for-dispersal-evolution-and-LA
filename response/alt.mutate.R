## mutation process

## A-locus
mutate.A <- function(prms, pop){
  
  ## rbinom mutants with prob=mr
  num.A.mutants <- rbinom(1, sum(pop), prob=prms$mr.A)
  ## break out if no mutants
  if(num.A.mutants==0) return(pop)
  mut.pop <- array(rmultinom(1, num.A.mutants, prob=pop),
                   dim=dim(pop))
  ## make sure not more mutants than individuals in each patch
  mut.pop[mut.pop>pop] <- pop[mut.pop>pop]
  
  ind.from <- which(mut.pop>0, arr.ind=TRUE)
  ## number of mutants per patch
  num.mutants <- mut.pop[ind.from]

  ## move from [,,1,] to [,,2,] or vice versa
  ind.to <- ind.from
  ind.to[,'dim3'] <- ind.from[,'dim3',drop=FALSE]%%2+1
  pop[ind.from] <- pop[ind.from]-num.mutants
  pop[ind.to] <- pop[ind.to]+num.mutants
  pop
}

## dispersal strategy
mutate.strategy.v2 <- function(prms, pop, strategies) {

  if(nrow(strategies)>=prms$num.strategies)
    return(list(pop=pop, strategies=strategies))

  max.num.mut <- min(prms$num.strategies-nrow(strategies), sum(pop))
  mutants <- array(rmultinom(1,max.num.mut,prob=pop),
                   dim=dim(pop))

  ## make sure no population has more mutants than individuals
  too.many <- mutants>pop
  mutants[too.many] <- pop[too.many]
  num.mut <- sum(mutants)

  ## check for possible problem
  if(any(mutants>pop)) cat('ERROR: too many mutants\n')

    ## remove mutants
  pop <- pop-mutants
  
  mutant.ind <- which(mutants>0, arr.ind=TRUE)
  mutant.ind <- mutant.ind[rep(1:nrow(mutant.ind),
                               mutants[mutants>0]),,drop=FALSE]
  colnames(mutant.ind) <- names(dim(pop))
  
  ## create new strategies based on existing strategies
  ##
  ## strategies that are going to mutate:
  mut.start <- strategies[mutant.ind[,'M'],, drop=FALSE]
  
  ## Pick 2 bins from each strategy
  mut.bins <- t(replicate(num.mut, sample(prms$n.bins+1, 2)))
  mut.loc <- cbind(1:num.mut, mut.bins)
  colnames(mut.loc) <- c('ind', 'b1', 'b2')

  ## effect sizes
  mut.sums <-
    mut.start[mut.loc[,c('ind','b1'),drop=FALSE]] +
    mut.start[mut.loc[,c('ind','b2'),drop=FALSE]]

  ## individuals with less than 0.01 mutational potential
  zero.cases <- which(mut.sums<=0.01)
  num.sum.to.zero <- length(zero.cases)
  
  if(num.sum.to.zero>0) {
    ## draw non-zero bin 2 for non-keepers; don't have to worry about
    ## accidentally drawing the same bin twice because we know bin 1
    ## is 0 in all cases.
    mut.start.zero.cases <- mut.start[zero.cases,, drop=FALSE]

    new.bins <- sapply(1:num.sum.to.zero,
                       function(x)
                         sample(which(mut.start.zero.cases[x,]>0.01), size=1))

    mut.loc[zero.cases,'b2'] <- new.bins
    
    ## effect sizes
    mut.sums <-
      mut.start[mut.loc[,c('ind','b1'),drop=FALSE]] +
      mut.start[mut.loc[,c('ind','b2'),drop=FALSE]]
  }

  ## enforce minimum effect size
  mut.effect <- runif(num.mut, min=0.1, max=0.9)*mut.sums

  ## redistribute
  new.strategies <- mut.start
  new.strategies[mut.loc[,c('ind','b1'),drop=FALSE]] <- mut.effect
  new.strategies[mut.loc[,c('ind','b2'),drop=FALSE]] <- mut.sums-mut.effect 
  strategies <- rbind(strategies, new.strategies)
   
  if(nrow(strategies) != nrow(unique(strategies)))
    cat('!!!!!!!!!!!!!!!!!!!!Duplicate strategies found!!!!!!!!!!!!!!!!!!!!\n')

  ## extend fourth dimension of pop array
  ##
  ## create mutants
  pop.mutants <- array(0, dim=c(dim(pop)[c('x', 'y', 'A')],
                                M=nrow(new.strategies)))
  new.indices <- mutant.ind
  new.indices[,4] <- 1:num.mut
  pop.mutants[new.indices] <- 1

  ## add mutants to original pop
  new.pop <- abind(pop, pop.mutants, along=4)

  names(dim(new.pop)) <- names(dim(pop))
  list(pop=new.pop, strategies=strategies)
}

## dispersal strategy
mutate.strategy <- function(prms, pop, strategies) {

  if(nrow(strategies)>=prms$num.strategies)
    return(list(pop=pop, strategies=strategies))

  max.num.mut <- min(prms$num.strategies-nrow(strategies), sum(pop))
  mutants <- array(rmultinom(1,max.num.mut,prob=pop),
                   dim=dim(pop))

  ## make sure no population has more mutants than individuals
  too.many <- mutants>pop
  mutants[too.many] <- pop[too.many]
  num.mut <- sum(mutants)

  ## check for possible problem
  if(any(mutants>pop)) cat('ERROR: too many mutants\n')
  
  ## remove mutants
  pop <- pop-mutants
  
  mutant.ind <- which(mutants>0, arr.ind=TRUE)
  mutant.ind <- mutant.ind[rep(1:nrow(mutant.ind),
                               mutants[mutants>0]),,drop=FALSE]
  colnames(mutant.ind) <- names(dim(pop))
  
  ## create new strategies based on existing strategies
  ##
  ## strategies that are going to mutate:
  mut.start <- strategies[mutant.ind[,'M'],, drop=FALSE]

  ## create dd (probability of staying home) and drop first column
  ## from mut.start 
  dd <- mut.start[,1,drop=FALSE]
  mut.start <- mut.start[,-1,drop=FALSE]
  
  ## mutate d and make sure values are between 0.001 and 0.999 (note -
  ## letting them equal 0.0001 or 0.999 can cause problems later
  ## on...)
  
  #MLC: alternate method of calculating mutational effect, enforcing min and max effect size:
  mut.effect<-runif(length(dd), min=0.5*prms$mut.size, max=(0.5*prms$mut.size+prms$mut.size))
  mut.direction<-sample(c(1,-1), length(dd), replace=TRUE, prob=c(0.5,0.5))
  dd <- dd+mut.effect*mut.direction
  
  #mutational effect drawn from a normal with sd = mut.size 
  #dd <- rnorm(length(dd), mean=dd, prms$mut.size)
  dd <- pmin(pmax(dd, 0.0001), 0.9999)

  ## calculate 'kk' for each strategy
  kk <- apply(mut.start!=0, 1, function(x) max(which(x)))

  ## mutate kk
  kk <- kk + sample(c(-1,0,1), size=length(kk), replace=TRUE)
  
  ## MLC: Being restricted to single step mutations may keep pops where long distance, but not short distance, dispersal is favored from crossing over the fitness valley of high frequency short distance dispersal to low frequency long distance dispersal.
  ## alternative kk mutation, with possibility to gain or lose up to max.dist/2:
#  dist.options<-seq(-0.5*prms$max.dist, 0.5*prms$max.dist)
#  probs<-c(seq(1,50,length=(length(dist.options)/2-1)), 50, seq(50,1, length=(length(dist.options)/2-1)))
#  kk<-kk + sample(dist.options, size=length(kk), replace=TRUE, prob=probs)

  kk <- pmin(pmax(kk, 1), ncol(mut.start))

  ## fill in values 1:kk[i] for each strategy, i
  mut.start <- mut.start*0
  for(i in 1:length(kk)){
    mut.start[i,1:kk[i]] <- (1-dd[i])/kk[i]
  }

  mut.start <- cbind(dd, mut.start)

  ## redistribute
  new.strategies <- mut.start
  strategies <- rbind(strategies, new.strategies)
  
  ## extend fourth dimension of pop array
  ##
  ## create mutants
  pop.mutants <- array(0, dim=c(dim(pop)[c('x', 'y', 'A')],
                                M=nrow(new.strategies)))
  new.indices <- mutant.ind
  new.indices[,4] <- 1:num.mut
  pop.mutants[new.indices] <- 1

  ## add mutants to original pop
  new.pop <- abind(pop, pop.mutants, along=4)

  names(dim(new.pop)) <- names(dim(pop))
  list(pop=new.pop, strategies=strategies)
}
