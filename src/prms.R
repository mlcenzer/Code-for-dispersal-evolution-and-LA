## ******************************************************************
## simulation parameter set
## ******************************************************************
##
## n.gens: number of iterations
## print.every: interval which output is printed to R console
## file.name: name of file to save output to

make.prms.sim <- function(n.gens=100,
                          save.freq=NA,
                          print.every=as.integer(1),
                          initial.strat='uniform',
                          file.name='res.RData') {
  if(is.na(save.freq)) save.freq <- rounf(n.gens/10)
  list(n.gens=n.gens,
       save.freq=save.freq,
       print.every=print.every,
       initial.strat=initial.strat,
       file.name=file.name)
}
## ******************************************************************

## ******************************************************************
## landscape parameter set
## ******************************************************************

make.prms.ls <- function(ltc.size=128, ## grid size
                         acl=c(0.1,0.1), ## landscape autocorrelation
                                         ## for each habitat
                         frac.suitable=c(0.5,0)) {
  
  ## make landscape
  landscape <- make.landscape(ltc.size=ltc.size,
                              acl=acl,
                              frac.suitable=frac.suitable)
  inhabitable <- apply(landscape, 1:2, sum)>0

  list(ltc.size=ltc.size,
       acl=acl,
       frac.suitable=frac.suitable,
       landscape=landscape,
       inhabitable=inhabitable)
}
## ******************************************************************

## ******************************************************************
## model parameter set
## ******************************************************************

make.prms <- function(landscape,
                      max.dist=10, ## max dispersal distance
                      n.bins=20, ## number of bins in dispersal kernel
                      num.strategies=50,
                      n.alleles=c(A=2, M=1), ## (A-locus, M-locus)
                      r=1.01,
                      K=10,
                      mut.size=0.01,
                      recomb=0,
                      ws=c(A.h1=1, a.h1=1, A.h2=1, a.h2=1),
                      mr.A=1e-4,
                      disturb.r=0,
                      ...) {

  if(n.bins>max.dist) {
    cat('ERROR: n.bins>max.dist\n')
    break
  }
  
  dn <- list(genotype=c('A', 'a'),
             habitat=c('h1','h2'))

  fitness <- matrix(ws,
                    nrow=n.alleles['A'],
                    ncol=2, ## REPLACE 2 WITH DIM OF LANDSCAPE
                    dimnames=dn)
  
  ## ------------------------------------------------------------
  ## create distance class object for possible dispersal steps
  
  step.dist.list <- make.step.dist(max.dist=max.dist, n.bins=n.bins)
  step.dist  <- step.dist.list[[1]]
  dist.class <- step.dist.list[[2]]
  
  ## ------------------------------------------------------------

  ## inputs <- lapply(as.list(match.call())[-1], eval.parent)
  ## inputs.f <- lapply(formals(), eval)
  ## for(v in names(inputs.f))
  ##   if(!(v %in% names(inputs)))
  ##     inputs <- append(inputs, inputs.f[v]) ## add formal args
  ## inputs[order(names(inputs))] ## order for consistency

  ## make fitness / growth / etc matrics
  ##
  ## because this matrix is the same for all three cases below
  ## (fitness, r, K), make it once and use it for all three
  null.mat <- array(0, dim=c(dim(landscape)[1],
                             dim(landscape)[2],
                             n.alleles['A'], 2))
  dimnames(null.mat) <- list(x=NULL,
                             y=NULL,
                             genotype=c('A', 'a'),
                             habitat=c('h1', 'h2'))

  ## create a function to populate null matrix with values from
  ## smaller matrix and landscape, rather than duplicating code to do
  ## it three times
  fill.matrix <- function(mat, f) {
    null.mat[,,'A','h1'][landscape[,,1]==1] <- mat['A','h1']
    null.mat[,,'a','h1'][landscape[,,1]==1] <- mat['a','h1']
    null.mat[,,'A','h2'][landscape[,,2]==1] <- mat['A','h2']
    null.mat[,,'a','h2'][landscape[,,2]==1] <- mat['a','h2']
    apply(null.mat, MARGIN=c(1,2,3), FUN=f)
  }

  ## Fitness for each genotype in each patch, assuming they use the
  ## habitat type in that patch with highest fitness; could also
  ## sample randomly.
  ##
  ## Survival is drawn as a binomial for each A-genotype from this
  ## array during selection.
  real.hab.fit <- fill.matrix(fitness, f=max)
  names(dim(real.hab.fit)) <- c('x', 'y', 'A')

  ## r by habitat & A genotype: currently r is the same everywhere &
  ## for each genotype, but we might want to change that later. Could
  ## be expanded to include r that varies by dispersal type.
  ##
  ## r for each genotype in each patch, assuming they use the habitat
  ## type in that patch with highest r; could also sample randomly, or
  ## maximize some other parameter (eg, survival)
  
  #real.hab.r <- fill.matrix(r.patch, f=max)

  ## K by habitat & A genotype: currently K is the same everywhere &
  ##for each genotype, but we might want to change that later.
  ##
  ## We should talk about how we want to implement K. If two habitat
  ## types both occur in a patch, does density regulation act
  ## independently on each type? (currently how real.hab.K is
  ## structured)
  #real.hab.K <- fill.matrix(K.patch, f=sum)

  ## whatever is needed from the below can go here
  ## real.hab.K<-apply(habitat.K, MARGIN=c(1,2,3), FUN=sum)
  ## names(dim(real.hab.K)) <- c('x', 'y', 'A')
  ## use this object later for decomposing into patches with one
  ## vs. two habitat types 
  ## shared.patches<-which(real.hab.K>max(K.patch))
  #real.hab.K[,,1][which(real.hab.K[,,genotype='A']>max(K.patch))] <-
  #  K.patch['A','h1']
  #real.hab.K[,,2][which(real.hab.K[,,genotype='a']>max(K.patch))] <-
  #  K.patch['a','h2']

  list(dim.pop=c(x=dim(landscape)[1], y=dim(landscape)[2], n.alleles),
       max.dist=max.dist,
       n.bins=n.bins,
       dist.class=dist.class,
       step.dist=as.matrix(step.dist),
       num.strategies=num.strategies,
       n.alleles=n.alleles,
       r=r,
       K=K,
       fitness=fitness,
       real.hab.fit=real.hab.fit,
       mut.size=mut.size,
       recomb=recomb,
       mr.A=mr.A,
       disturb.r=disturb.r)
}
## ******************************************************************
