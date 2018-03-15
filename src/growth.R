## birth: create offspring to found next generation
growth <- function(prms, pop, patch.sizes) {
  patch.sizes <- apply(pop, 1:2, sum)
  growth.rate <- array(1+prms$r*(1-patch.sizes/prms$K), dim=dim(pop))
  growth.rate[growth.rate<0] <- 0
  ## poisson number of offspring
  num.offspring <- rpois(n=length(pop), lambda=pop * growth.rate)
  offspring <- array(num.offspring, dim=dim(pop))
  offspring
}




## VERSION TO BE EXPANDED LATER
## patch.sizes <- apply(pop, c(1,2,3), sum)

## ## calculate expected number of offspring per type per patch, using
## ## a logistic model with growth rate r.patch and carrying capacity
## ## K.patch.  Currently, intra- and inter-type competition are
## ## equally strong.
## mean.num.offspring <-
##   patch.sizes*(1+prms$real.hab.r*(1-(patch.sizes/prms$real.hab.K)))

## ## LKM: added below line to deal with patches with K=0
## mean.num.offspring[is.nan(mean.num.offspring)] <- 0
## mean.num.offspring[mean.num.offspring<0] <- 0
## ## number of offspring per type per patch
## offspring <- apply(mean.num.offspring, MARGIN=c(1,2,3),
##                    FUN=function(x) rpois(1,x))

## ## LKM: below is my vector based version - its almost the same, but
## ## for some reason the apply statement above is dropping the M
## ## component of the population.  Is that what we want?
## offspring <- array(rpois(length(pop),mean.num.offspring),
##                    dim=dim(pop))

## ## this returns warnings; they aren't actually a problem (empty
## ## patches return an NA). This is probably avoidable, but if we make
## ## empty patches 0 in the mean.num.offspring array we might get
## ## spontaneous generation.
## ##
## ## LKM: I don't think we'd get spontaneous generation, as a poisson
## ## with mean 0 will always give zero.
## offspring[which(is.na(offspring))]<-0
## offspring
## ## needs to include dispersal genotype structure
