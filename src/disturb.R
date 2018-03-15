disturb <- function(prms, pop) {
  ## sample patches to 'disturb'
  disturbed <- array(rbinom(length(pop), 1, prob=1-prms$disturb.r),
                     dim=dim(pop))
  pop*disturbed
}
