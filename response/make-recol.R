## ***************************************************
setwd('~/Dropbox/sbb_dispersal/displa')
rm(list=ls())
source('src/initialize.R')
setwd('~/Documents/sbb_dispersal/displa') #remove for files in Dropbox
## ***************************************************

make.summary <- function(save.path, dir, n.cores=1) {
  load(file.path(save.path, 'prms.array.RData'))
  load(file.path(save.path, 'prms.ls.list.RData'))
  completed <- list.files(file.path(save.path, dir))
	
  get.recol.vectors <- function(x) {
    cat(x, '\n')
    load(file.path(save.path, dir, x))
    if(is.null(out)) return(NA)
    pop.list   <- out$pop.list
    strat.list <- out$strat.list
	recol.list <- out$recol.list

    run <- as.numeric(gsub('.RData', '', x))
    H0 <- prms.ls.list[[run]]$landscape
    H <- abind(H0, H0[,,1]*H0[,,2], along=3) #all H1, all H2, patches with both
    H[,,1] <- H0[,,1] - H[,,3]
    H[,,2] <- H0[,,2] - H[,,3] 
    
		recol.array<-array(unlist(recol.list), dim=c(dim(recol.list[[1]]),length(recol.list)))[,,2:length(recol.list)]

		#convert recol states into strings through time for each patch
		recol.strings<-apply(recol.array, 1:2, function(x) paste(x, collapse=''))
		
		#split recol.strings by patch type here
		H1.recol.strings<-recol.strings
			H1.recol.strings[which(H[,,1]==0)]<-''
		H2.recol.strings<-recol.strings
			H2.recol.strings[which(H[,,2]==0)]<-''
		H3.recol.strings<-recol.strings
			H3.recol.strings[which(H[,,3]==0)]<-''
			
		#############################		
		recol.by.string<-function(string){
			times<-c()
			events.list<-strsplit(paste(string, 'R'), 'R') #save cases that end with recol
			events<-unlist(events.list)
			if(length(events)==1){
				return(times) #return empty vector; no recolonizations
				}
			#drop last element
			events<-events[1:(length(events)-1)]			
			for(n in 1:length(events)){
				times[n]<-1+sum(as.numeric(unlist(strsplit(events[n], ''))))
				}
		times
		}
		
		recol.times.H1<-unlist(sapply(H1.recol.strings, recol.by.string), use.names=FALSE)
		recol.times.H2<-unlist(sapply(H2.recol.strings, recol.by.string), use.names=FALSE)
		recol.times.H3<-unlist(sapply(H3.recol.strings, recol.by.string), use.names=FALSE)
	
		recol.times<-list(recol.times.H1, recol.times.H2, recol.times.H3)
      
	list(recol.times=recol.times)
  }

  if(n.cores==1)
    recol.list <- lapply(completed, get.recol.vectors)
  if(n.cores>1)
    recol.list <- mclapply(completed, get.recol.vectors, mc.cores=n.cores)
                             
    names(recol.list) <- completed

  save(recol.list,
       file=sprintf('~/Dropbox/sbb_dispersal/displa/%s/%s.RData', save.path, dir))
}           
                                                                                                   

make.summary(save.path='saved/long runs/h1.acl-0.01;h1.frac.s-0.3',
             dir='disturb.r-0.1;recomb-0.0;initial.strat-5;s-0.50.a.and.A_recol_flatstrat', n.cores=1)

