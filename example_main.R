#####Set working directory, clear history, load functions
setwd('~/Dropbox/sbb_dispersal/displa')
rm(list=ls())
source('src/initialize.R')

#Example for generating landscapes with different parameter ranges
var.acl.frac.suitable.h2.only(acl.h1=0.01, #spatial autocorrelation in resource 1
                              frac.suitable.h1=0.3, #abundance of resource 1 (as fraction of landscape)
                              range.acl.h2=c(0.001,0.01,0.1), #range of spatial autocorrelations for resource 2
                              range.frac.suitable.h2=c(0.1,0.3,0,5) #range of abundances for resource 2
                              n.reps=25, #Number of replicate landscapes of each type
                              save.path='saved/response/init.strat', #where to save landscapes
                              remove.files=TRUE, #remove old files in save location
                              n.cores=25, #number of cores to run in parallel
                              ltc.size=128) #size of landscape grid
                              
run.prms.ls.list(save.path='saved/response/init.strat', #where to save simulation output
save.dir='disturb.r-0.1;recomb-0.0;initial.strat-5;s-0.50.a.and.A', #what to name simulation output
                 initial.strat=5, #starting strategy for the population
                 remove.files=TRUE, #remove old files
                 n.gens=2e4, #number of generations to run simulation
                 save.freq=2e2, #how often to save population
                 n.cores=25, #number of cores to run in parallel
                 n.bins=10, #number of bins in the dispersal kernel
                 max.dist=10, #maximum distance (note that the size of each bin = max.dist/n.bins)
                 disturb.r=0.1, #probability of patch extinction
                 ws=c(A.h1=1, #survival of ecological fit allele A in resource 1
                      a.h1=0.5, #survival of ecological fit allele a in resource 1
                      A.h2=0.5, #survival of ecological fit allele A in resource 2
                      a.h2=1), #survival of ecological fit allele a in resource 2
                 recomb=0) #recombination rate
                 #Note that other parameters may be added; see meta.functions for full range.