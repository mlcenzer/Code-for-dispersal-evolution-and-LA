next.gen <- function(prms, prms.ls, pop, strategies) {

  res <- list(pop=pop,
              strategies=strategies,
              md=NA,
              ms=NA,
              dg=NA)

  ## patch disturbance events
  pop.e <- disturb(prms, pop)
  if(sum(pop.e)==0) {
    cat('Extinction after disturbance.\n')
    res$pop <- pop.e
    return(res)
  }

  ## dispersal
  pop.d <- disperse(prms, prms.ls, pop.e, strategies)
  pop.d.2 <- dispersal.mortality(prms.ls, pop.d)
  pop.d.by.A <- apply(pop.d, 3, sum)
  pop.d.2.by.A <- apply(pop.d.2, 3, sum)
  res$md <- (pop.d.by.A-pop.d.2.by.A)/pop.d.by.A ## dispersal mortality
    if(sum(pop.d.2)==0) {
    cat('Extinction after dispersal.\n')
    res$pop <- pop.d.2
    return(res)
  }
  
  ## selection
  pop.s <- select(prms, pop.d.2)
  pop.s.by.A <- apply(pop.s, 3, sum)
  res$ms <- (pop.d.2.by.A-pop.s.by.A)/(pop.d.2.by.A) ## selection mortality
  if(sum(pop.s)==0) {
    cat('Extinction after selection.\n')
    res$pop <- pop.s
    return(res)
  }

  ## drop extinct strategies (good to do this before mutation)
  tmp <- loss(pop.s, strategies)
  pop.l <- tmp$pop
  strategies <- tmp$strategies
  
  ## mutation
  pop.m <- mutate.strategy(prms, pop.l, strategies)
  strategies <- pop.m$strategies  
  pop.m <- pop.m$pop
  pop.m.2 <- mutate.A(prms=prms, pop=pop.m)

  ## recombination
  pop.r <- recombine(prms=prms, pop=pop.m.2)
  pop.r.by.A <- apply(pop.r, 3, sum)
 
  ## growth
  pop.g <- growth(prms=prms, pop=pop.r)
  pop.g.by.A <- apply(pop.g, 3, sum)
  res$dg <- (pop.r.by.A-pop.g.by.A)/pop.r.by.A ## growth mortality
  if(sum(pop.g)==0) {
    cat('Extinction after growth.\n')
    res$pop <- pop.g
    res$strategies <- strategies
    return(res)
  }

  res$pop <- pop.g
  res$strategies <- strategies
  res
}
