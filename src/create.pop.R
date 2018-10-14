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

#####To add: is there a way to do this from a summary file?

load.pop <- function(ii, load.dir, save.path, load.strat, ...){

	load(file=sprintf('%s/%s/%s.RData', save.path, load.dir, ii))

	prev.pop<-out$pop.list[[length(out$pop.list)]]

	if(load.strat==FALSE) {nom.pop<-apply(prev.pop, 1:3, sum)
	return(array(nom.pop, dim=prms$dim.pop))}
	
	else return(prev.pop)
}

load.strategies <- function(ii, load.dir, save.path, initial.strat, load.strat, ...){
	load(file=sprintf('%s/%s/%s.RData', save.path, load.dir, ii))
	final_strat <- out$strat.list[[length(out$strat.list)]]
	if(load.strat==FALSE){
		flat_strat<-create.strategies(prms, initial.strat)
		return(flat_strat)
	}
	else return(final_strat)
}
