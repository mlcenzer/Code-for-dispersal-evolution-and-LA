loss <- function(pop, strategies ) {
  
  ## check for loss of a genotype
  lost <- which(apply(pop, 4, sum)==0)
  if(length(lost)==0) return(list(pop=pop, strategies=strategies))

  list(pop=pop[,,,-lost,drop=FALSE],
       strategies=strategies[-lost,,drop=FALSE])
}
