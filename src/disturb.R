disturb <- function(prms, pop) {
  ## sample patches to 'disturb'
  
  ####OH NO FIX ME!
  #disturbed <- array(rbinom(length(pop), 1, prob=1-prms$disturb.r),
  #                   dim=dim(pop))
  
  disturbed <- array(rbinom(length(pop[,,1,1]), 1, prob=1-prms$disturb.r),
                     dim=dim(pop))
  
 pop*disturbed
}
