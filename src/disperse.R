## individuals disperse
disperse <- function(prms, prms.ls, pop, strategies) {
  
  ## first, figure out who disperses
  prob.leave <- rep(1-strategies[,1], each=prod(dim(pop)[c('x','y','A')]))
  disp <- array(rbinom(length(pop), size=pop, prob=prob.leave),
                dim=dim(pop), dimnames=dimnames(pop))
                
  if(sum(disp)==0) return(pop=pop)
  
  pop <- pop-disp # remove the dispersers from pop 

  ## second, figure out where dispersers go

  ## create array of disperser indices
  tmp <- which(disp>0, arr.ind=TRUE)
  disp.ind <- tmp[rep(1:nrow(tmp), disp[disp>0]),,drop=FALSE]
  colnames(disp.ind) <- names(dim(pop))
  ## drop first column of strategies because we are now only focussing
  ## on dispersers; index by M-genotype to extract correct dispersal
  ## function for each disperser
  disp.strat <- strategies[disp.ind[,'M'],-1,drop=FALSE]

  ## generate dispersal distances
  disp.ind <- cbind(disp.ind,
                    dest.bin=rep(0, nrow(disp.ind)),
                    dest.dist.class.row=rep(0, nrow(disp.ind)))
  for(strategy in unique(disp.ind[,'M'])) {
    inds <- disp.ind[,'M']==strategy
    disp.ind[inds,'dest.bin'] <- 
      sample.int(prms$n.bins,
                 size=sum(inds),
                 prob=strategies[strategy,-1],
                 replace=TRUE)
  }
  for(dest.bin in unique(disp.ind[,'dest.bin'])) {
    inds <- disp.ind[,'dest.bin']==dest.bin
    disp.ind[inds,'dest.dist.class.row'] <- 
      sample(prms$dist.class[[dest.bin]],
             size=sum(inds),
             replace=TRUE)
  }
  disp.rows <- disp.ind[,'dest.dist.class.row']
  
  ## determine steps in the x and y directions and then update
  ## disp.ind accordingly
  xy <- c('x','y')
  xyAM <- c('x','y','A','M')
  delta.xy <- prms$step.dist[disp.rows, xy]
  disp.ind[,xy] <- (disp.ind[,xy] + delta.xy) %% prms.ls$ltc.size
  disp.ind[,xy][disp.ind[,xy]==0] <- prms.ls$ltc.size

  ## add dispersers in new spots - this should be looked at again for
  ## speeding up - if we only had one of each type of each genotype
  ## arriving in each patch it would be much faster.
  ##
  for(i in 1:nrow(disp.ind))
   pop[disp.ind[i,xyAM,drop=FALSE]] <- pop[disp.ind[i,xyAM,drop=FALSE]]+1
  
  pop
}

dispersal.mortality <- function(prms.ls, pop) {
  inhabitable <- array(prms.ls$inhabitable, dim=dim(pop))
  pop*inhabitable
}
