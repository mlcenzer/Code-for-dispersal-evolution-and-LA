## ***************************************************
setwd('~/Dropbox/sbb_dispersal/displa')
rm(list=ls())
source('src/initialize.R')
#setwd('~/Documents/sbb_dispersal/displa') #remove for files in Dropbox
## ***************************************************

make.summary <- function(save.path, dir, n.cores=1) {
  load(file.path(save.path, 'prms.array.RData'))
  load(file.path(save.path, 'prms.ls.list.RData'))
  completed <- list.files(file.path(save.path, dir))

  get.summary.vectors <- function(x) {
    cat(x, '\n')
    load(file.path(save.path, dir, x))
    if(is.null(out)) return(NA)
    pop.list   <- out$pop.list
    strat.list <- out$strat.list

    get.geno.freqs <- function(pop) {
      if(sum(pop)==0 | is.null(pop)) return(NA)
      apply(pop, 3:4, sum)
    }
    get.geno.space.freqs <- function(pop){
      if(sum(pop)==0 | is.null(pop)) return(NA)
      apply(pop, 1:3, sum)    	
    }
    get.a.hab.freqs <- function(pop, landscape, landscape.sums){
      #browser()
      if(is.null(pop)) return(NA)
      else if(is.na(pop)) return(NA)
      else if(length(dim(pop))!=3) browser()
      else{ 
      	A1<-abind(pop[,,1], pop[,,1], pop[,,1], along=3)
        A2<-abind(pop[,,2], pop[,,2], pop[,,2], along=3)
        HA<-A1*H
        Ha<-A2*H
        full.freqs <- abind(HA, Ha, along=4)
        dimnames(full.freqs)[[3]] <- c("H1", "H2", "H3") 
    	dimnames(full.freqs)[[4]] <- c("A", "a")
    	Asums<-apply(full.freqs, 3:4, sum)
    	Asums/landscape.sums}
    }
    
    geno.freqs <- lapply(pop.list, get.geno.freqs)
    a.freqs <- lapply(geno.freqs, function(x)
      if(any(is.na(x))) { x } else { apply(x, 1, sum) })
    m.freqs <- lapply(geno.freqs, function(x)
      if(any(is.na(x))) { x } else { apply(x, 2, sum) })

                                        #get A distribution
    run <- as.numeric(gsub('.RData', '', x))
    H0 <- prms.ls.list[[run]]$landscape
    H <- abind(H0, H0[,,1]*H0[,,2], along=3) #all H1, all H2, patches with both
    H[,,1] <- H0[,,1] - H[,,3]
    H[,,2] <- H0[,,2] - H[,,3]    
    Hsums<-apply(H, 3, sum)
	#if(run==37) browser()
    a.space.freqs <- lapply(pop.list, get.geno.space.freqs)
    a.freqs.space.time <- lapply(a.space.freqs, get.a.hab.freqs, landscape=H, landscape.sums=Hsums)

    ## patch.sharing<-lapply(a.space.freqs, function(x){
    ##   x[which(x>0)]<-1
    ## 					#x[,,1]+x[,,2]
    ##   x
    ## })

    make.m.by.h<-function(pop, strat, h) {
      if(length(pop)==0) return(rep(NA, prms$n.bins+1))
      if(is.null(pop)) return(NA)
      if(is.na(pop)) return(NA)
      habs3<-prms.ls.list[[run]]$landscape[,,1]*prms.ls.list[[run]]$landscape[,,2]
      habs<-abind(prms.ls.list[[run]]$landscape[,,1]-habs3,
                  prms.ls.list[[run]]$landscape[,,2]-habs3, habs3, along=3)
      
      H<-array(rep(habs[,,h], times=(dim(pop)[3])*(dim(pop)[4])), dim=dim(pop))
      m<-apply(pop*H, 4, sum)
      m.h<-colSums(m*strat)/sum(m)
     m.h
    }

    m.by.h1<-sapply(1:length(pop.list), FUN= function(x) make.m.by.h (pop=pop.list[[x]], strat=strat.list[[x]], h=1))
    m.by.h2<-sapply(1:length(pop.list), FUN= function(x) make.m.by.h (pop=pop.list[[x]], strat=strat.list[[x]], h=2))
    m.by.h3<-sapply(1:length(pop.list), FUN= function(x) make.m.by.h (pop=pop.list[[x]], strat=strat.list[[x]], h=3))
    m.by.h<-array(c(m.by.h1,m.by.h2, m.by.h3), dim=c(dim(strat.list[[1]])[2],length(pop.list), 3))
    
    list(geno.freqs=geno.freqs,
         a.freqs=a.freqs,
         m.freqs=m.freqs,
         a.space.freqs=a.freqs.space.time,
                                        #         patch.sharing=patch.sharing,
         strats=strat.list,
         mds.vec=out$mds,
         mss.vec=out$mss,
         deltags.vec=out$deltags,
         m.by.h=m.by.h)
  }

  if(n.cores==1)
    summary.list <- lapply(completed, get.summary.vectors)
  if(n.cores>1)
    summary.list <- mclapply(completed, get.summary.vectors,
                             mc.cores=n.cores)

  make.list <- function(s) {
    tmp <- lapply(summary.list, function(x) x[[s]])
    names(tmp) <- completed
    tmp
  }

  geno.freqs.list <- make.list('geno.freqs')
  a.freqs.list <- make.list('a.freqs')
  m.freqs.list <- make.list('m.freqs')
  a.space.freqs.list <- make.list('a.space.freqs')
                                        #  patch.sharing <- make.list('patch.sharing')
  strats.list <- make.list('strats')
  mds.list <- make.list('mds.vec')
  mss.list <- make.list('mss.vec')
  deltags.list <- make.list('deltags.vec')
  m.by.h.list <- make.list('m.by.h')

  save(geno.freqs.list,
       a.freqs.list,
       m.freqs.list,
       a.space.freqs.list,
                                        #       patch.sharing,
       strats.list,
       mds.list,
       mss.list,
       deltags.list,
       m.by.h.list,
       file=sprintf('~/Dropbox/sbb_dispersal/displa/%s/%s.RData', save.path, dir))
}           

#make.summary(save.path='saved/response/init.strat', dir='disturb.r-0.1;recomb-0.0;initial.strat-5;s-0.50.a.and.A;init.strat-1', n.cores=75)#I am so sorry this is what this is called.
 
                                                                    
make.summary(save.path='saved/response/K', dir='disturb.r-0.1;recomb-0.0;initial.strat-5;s-0.50.a.and.A.K-100', n.cores=75)

