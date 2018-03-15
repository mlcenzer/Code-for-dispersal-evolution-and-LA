create.pop <- function(prms, prms.ls) {
  ## begin with every habitable cell occupied at K
  pop.sum <- (prms.ls$landscape[,,1] | prms.ls$landscape[,,2])*prms$K
  ## works for 2 A alleles & 1 starting M
  pop.A <- apply(pop.sum, 1:2, function(y) sample(x=0:y, size=1))
  pop.a <- pop.sum - pop.A
  popAs <- abind(pop.A, pop.a, along=3)
  array(popAs, dim=prms$dim.pop)
}

create.strategies <- function(prms, initial.strat=1) {
  if(initial.strat>prms$n.bins) initial.strat <- prms$n.bins
  strategies <- matrix(0, nrow=prms$dim.pop['M'], ncol=prms$n.bins+1)
  strategies[1:initial.strat] <- 1/initial.strat
  strategies
}
