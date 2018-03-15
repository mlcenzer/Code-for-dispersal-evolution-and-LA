## perform recombination
recombine <- function(prms, pop){

  ## no recombination if recomb rate is 0
  if(prms$recomb==0) return(pop)
  ## no recombination if only one allele at A-locus
  if(any(apply(pop, 3, sum)==0)) return(pop)

  patch.sizes <- apply(pop, 1:2, sum)
  num.sex <- rbinom(n=length(patch.sizes),
                    size=floor(patch.sizes/2),
                    prob=prms$recomb)
  dim(num.sex) <- dim(patch.sizes)
  if(all(num.sex==0)) return(pop)

  recomb.patches <- which(num.sex>0, arr.ind=TRUE)
  colnames(recomb.patches) <- c('x', 'y')

  ## possible genotypes (could go in prms)
  sample.ind <- as.matrix(expand.grid(1:dim(pop)['A'], 1:dim(pop)['M']))

  get.parents <- function(x, y) {
    num <- num.sex[x,y]
    parents <- sample(rep(1:nrow(sample.ind), pop[x,y,,]),
                      size=2*num)
    list(mums=parents[1:num], dads=parents[(1:num)+num])
  }
  parents.ind <- mapply(get.parents,
                        x=recomb.patches[,1,drop=FALSE],
                        y=recomb.patches[,2,drop=FALSE],
                        SIMPLIFY=FALSE)

  mums <- lapply(parents.ind, function(x) x[['mums']])
  dads <- lapply(parents.ind, function(x) x[['dads']])
  mums.genotype <- sample.ind[unlist(mums),,drop=FALSE]
  dads.genotype <- sample.ind[unlist(dads),,drop=FALSE]
  colnames(mums.genotype) <- c('mum.A', 'mum.M')
  colnames(dads.genotype) <- c('dad.A', 'dad.M')

  ## create table (rows are matings)
  mate.table <- cbind(recomb.patches[rep(1:nrow(recomb.patches),
                                         num.sex[recomb.patches]),,
                                     drop=FALSE],
                      mums.genotype,
                      dads.genotype)

  ## drop rows where parents share an allele at either locus
  diff.at.both <-
    mate.table[,'mum.A']!=mate.table[,'dad.A'] &
    mate.table[,'mum.M']!=mate.table[,'dad.M']
  mate.table <- mate.table[diff.at.both,,drop=FALSE]
  if(nrow(mate.table)==0) return(pop)
  
  p1 <- c('x','y','mum.A','mum.M')
  p2 <- c('x','y','dad.A','dad.M')
  r1 <- c('x','y','mum.A','dad.M')
  r2 <- c('x','y','dad.A','mum.M')

  ## decrement parental genotypes
  pop[mate.table[,p1,drop=FALSE]] <- pop[mate.table[,p1,drop=FALSE]]-1
  pop[mate.table[,p2,drop=FALSE]] <- pop[mate.table[,p2,drop=FALSE]]-1
  ## increment recombinant genotypes
  pop[mate.table[,r1,drop=FALSE]] <- pop[mate.table[,r1,drop=FALSE]]+1
  pop[mate.table[,r2,drop=FALSE]] <- pop[mate.table[,r2,drop=FALSE]]+1

  pop
}

