run.sim <- function(prms,
                    prms.ls,
                    prms.sim,
                    print=FALSE,
                    plot=FALSE,
                    load.dir=NA,
                    ii,
                    save.path,
                    ...) {

#load pop from previous run
  if(!is.na(load.dir)){
  	pop <- load.pop(ii, load.dir, save.path, ...) 
	strategies <- load.strategies(ii, load.dir, save.path, initial.strat=prms.sim$initial.strat, ...)
	}
#create new pop
  else	{
  	pop <- create.pop(prms, prms.ls)
  	strategies <- create.strategies(prms, initial.strat=prms.sim$initial.strat)
  }

  recol <- array('c', dim=dim(pop)[1:2])

  n.saves <- ceiling(prms.sim$n.gens/prms.sim$save.freq)+1
  pop.list   <- vector('list', n.saves)
  strat.list <- vector('list', n.saves)
  recol.list <- vector('list', n.saves)
  mds <- matrix(nrow=2, ncol=n.saves)
  mss <- matrix(nrow=2, ncol=n.saves)
  deltags <- matrix(nrow=2, ncol=n.saves)
  md <- c(0,0) # this needs an initial value to get through the first if loop
  ms <- c(0,0)
  dg <- c(0,0)

  for(i in 0:prms.sim$n.gens) {

    if(i%%prms.sim$save.freq==0) {
      pop.list[[(i/prms.sim$save.freq)+1]]   <- pop
      strat.list[[(i/prms.sim$save.freq)+1]] <- strategies
      recol.list[[(i/prms.sim$save.freq)+1]] <- recol
      mds[,(i/prms.sim$save.freq)+1]     <- md
      mss[,(i/prms.sim$save.freq)+1]     <- ms
      deltags[,(i/prms.sim$save.freq)+1] <- dg
   }
    
    if(sum(pop)==0) {
      if(print) cat(sprintf('Population extinct at generation %d\n', i))
      return(list(prms=prms,
                  prms.ls=prms.ls,
                  prms.sim=prms.sim,
                  pop.list=pop.list,
                  strat.list=strat.list,
                  recol.list=recol.list,
                  extinction.time=i,
                  mds=mds,
                  mss=mss,
                  deltags=deltags))
    }
    
    ptm.1 <- proc.time()
    ng <- next.gen(prms, prms.ls, pop, strategies)
    ptm.2 <- proc.time()
    if(print) {
      cat(sprintf('%d: %2.3f seconds\n', i, (ptm.2 - ptm.1)[['elapsed']]))
      cat(sprintf('%d strategies, %d individuals\n',
                  nrow(ng$strategies),
                  sum(ng$pop)))
    }

    if(plot) {
      habitat <- cbind(rbind(which(prms.ls$landscape[,,1]==1,
                                   arr.ind=TRUE),
                             which(prms.ls$landscape[,,2]==1,
                                   arr.ind=TRUE)), 
                       hab=rep(c(1,2),
                               times=c(nrow(which(prms.ls$landscape[,,1]==1,
                                                  arr.ind=TRUE)),
                                       nrow(which(prms.ls$landscape[,,2]==1,
                                                  arr.ind=TRUE))))) 
      par(mfrow=c(1,2))
      plot(habitat[,1], habitat[,2],
           col=c(rgb(1,0,0,.5),rgb(0,0,1,.5))[habitat[,3]], pch=16)
      image(max(pop)-apply(pop, 1:2, sum), col=gray(level=0:10/10))
    }  
    pop <- ng$pop
    strategies <- ng$strategies
    recol <- ng$recol
    md <- ng$md
    ms <- ng$ms
    dg <- ng$dg
  } 
  
  return(list(prms=prms,
              prms.ls=prms.ls,
              prms.sim=prms.sim,
              pop.list=pop.list,
              strat.list=strat.list,
              recol.list=recol.list,
              extinction.time=NA,
              mds=mds,
              mss=mss,
              deltags=deltags))
}


var.acl.frac.suitable <- function(range.acl,
                                  range.frac.suitable,
                                  n.reps,
                                  save.path,
                                  remove.files,
                                  n.cores=1) {

  ## remove files from previous run
  if(remove.files==TRUE) {
    unlink(file.path(save.path), recursive=TRUE)
    dir.create(save.path)
  }

  fn.prms.array <- file.path(save.path, 'prms.array.RData')
  prms.array <- expand.grid(acl=range.acl,
                            frac.suitable=range.frac.suitable,
                            rep=1:n.reps)
  save(prms.array, file=fn.prms.array)

  make.prms.ii <- function(ii)
    prms.ls <- make.prms.ls(acl=c(prms.array$acl[ii], prms.array$acl[ii]),
                            frac.suitable=c(prms.array$frac.suitable[ii], 0))
  reps <- 1:nrow(prms.array)

  if(n.cores==1) prms.ls.list <- lapply(reps, make.prms.ii)
  if(n.cores>1) prms.ls.list <- mclapply(reps, make.prms.ii,
                                         mc.preschedule=FALSE,
                                         mc.cores=n.cores)
  save(prms.ls.list, file=file.path(save.path, 'prms.ls.list.RData'))
  NULL
}

var.acl.frac.suitable.h2.only <- function(acl.h1,
                                          frac.suitable.h1,
                                          range.acl.h2,
                                          range.frac.suitable.h2,
                                          n.reps,
                                          save.path,
                                          remove.files,
                                          n.cores=1,
                                          ltc.size=128) {

  ## remove files from previous run
  if(remove.files==TRUE) {
    unlink(file.path(save.path), recursive=TRUE)
    dir.create(save.path)
  }

  fn.prms.array <- file.path(save.path, 'prms.array.RData')
  prms.array <- expand.grid(acl=range.acl.h2,
                            frac.suitable=range.frac.suitable.h2,
                            rep=1:n.reps)
  save(prms.array, file=fn.prms.array)

  make.prms.ii <- function(ii, ltc)
    prms.ls <- make.prms.ls(ltc.size=ltc, acl=c(acl.h1,
                                  prms.array$acl[ii]),
                            frac.suitable=c(frac.suitable.h1,
                                            prms.array$frac.suitable[ii]))
  reps <- 1:nrow(prms.array)

  if(n.cores==1) prms.ls.list <- lapply(reps, make.prms.ii, ltc=ltc.size)
  if(n.cores>1) prms.ls.list <- mclapply(reps, make.prms.ii,
                                         mc.preschedule=FALSE,
                                         mc.cores=n.cores, ltc=ltc.size)
  save(prms.ls.list, file=file.path(save.path, 'prms.ls.list.RData'))
  NULL
}
  
run.prms.ls.list <- function(save.path,
                             save.dir,
                             initial.strat,
                             remove.files,
                             n.gens,
                             save.freq,
                             n.cores=1,
                             K=10,
                             num.strategies=50,
                             ...) {

## remove files from previous run
  if(remove.files==TRUE) {
    unlink(file.path(save.path, save.dir), recursive=TRUE)
    dir.create(file.path(save.path, save.dir))
  }
  
  ## load prms.ls.list
  load(file.path(save.path, 'prms.ls.list.RData'))
  nrep <- length(prms.ls.list)
  
  run.rep <- function(ii) {
    time.start <- Sys.time()
    prms.sim <- make.prms.sim(n.gens=n.gens,
                              save.freq=save.freq,
                              initial.strat=initial.strat)
    prms.ls <- prms.ls.list[[ii]]

## NORMAL make.prms
#    prms <- make.prms(landscape=prms.ls$landscape, ...)
## make.prms for review response with K=100
#	prms <- make.prms(landscape=prms.ls$landscape, K=100, ...)
## make.prms for review response with num.strategies=10 or 100 or whatever you please
	prms <- make.prms(landscape=prms.ls$landscape, num.strategies=num.strategies, K=K, ...)
	out <- run.sim(prms,
                   prms.ls,
                   prms.sim,
                   ii=ii,
                   save.path=save.path,
                   ...)
                   
    
    ## could extract a few things here if we don't want to keep everything
    save(prms, prms.ls, prms.sim, out,
         file=file.path(save.path, save.dir, sprintf('%d.RData', ii)))
    time.stop <- Sys.time()
    
    ## check how many replicates are complete
    files <- list.files(file.path(save.path, save.dir))
    cat(sprintf('%d/%d complete.  Replicate %d took %d minutes. ',
                length(files), nrep,
                ii, round(as.numeric(difftime(time.stop, time.start,
                                              units='mins')))))
    if(!is.na(out$extinction.time)) {
      cat(sprintf('Extinction at iteration %d.', out$extinction.time))
    }
    cat('\n')
  }
  
  completed <- list.files(file.path(save.path, save.dir))
  
  while(length(completed)<nrep) {
    rep.done <- sapply(strsplit(completed, '\\.'),
                       function(x) as.numeric(x[[1]]))
    to.do <- (1:nrep)[!(1:nrep %in% rep.done)]
    cat(sprintf('%d/%d remaining.\n', length(to.do), nrep))

    if(length(to.do)>n.cores & length(to.do)>1)
      to.do <- sample(to.do, replace=FALSE)
    if(n.cores==1) lapply(to.do, run.rep)
    if(n.cores>1) mclapply(to.do, run.rep,
                           mc.preschedule=FALSE,
                           mc.cores=n.cores)
    completed <- list.files(file.path(save.path, save.dir))
  }
  NULL
}
