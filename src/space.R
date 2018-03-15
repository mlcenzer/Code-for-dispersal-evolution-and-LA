## Habitat quality: array with 1,2 space and 3 the survival
## probability for each local adaptation genotype. 

## If we at some point add mutation in local adaptation genotype, we
## will want to change this structure to expand/contract as new
## genotypes arise. 

## landscape functions
locations.clumped <- function(size, phi, n.patch) {
  
  PosSymSqRt <- function(cov.mat, tol=10^-10) {
    eigs <- eigen ( cov.mat )
    keep <- eigs$values > tol
    rot.mat <- eigs$vectors[,keep] %*%
      diag(sqrt(eigs$values[keep])) %*%
      t(eigs$vectors[,keep]) 
    rot.mat
  }
  
  rmvn <- function(n, mu, V) {
    p <- length(mu)
    if(any(is.na(match(dim(V), p)))) 
      stop("Dimension problem!") 
    D <- PosSymSqRt(V, 10^-10) 
    t(matrix(rnorm(n * p), ncol = p) %*% D +
      rep(mu, rep(n, p)))
  }   

  dist.wa <- function(x) {
    f <- function(i, j)
      sqrt(min(abs(x[i,1] - x[j,1]), size - abs(x[i,1] - x[j,1]))^2 +
           min(abs(x[i,2] - x[j,2]), size - abs(x[i,2] - x[j,2]))^2)
    outer(1:nrow(x), 1:nrow(x), FUN=Vectorize(f))
  }
  
  ##code for evaluating clumpiness of underlying probability landscape      
  ##image(matrix(X, nrow=sqrt(nrow(simgrid))))
  ##library(ape)
  ##inv.dist<-1/clump.dist
  ##diag(inv.dist)<-0
  ##Moran.I(as.vector(X), inv.dist)

  simgrid <- expand.grid(seq(1, size, by=1), seq(1, size, by=1))
  clump.dist <- as.matrix(dist.wa(simgrid))
  X <- rmvn(1, rep(0, nrow(simgrid)), exp(-phi * clump.dist))
  weighted.pts <- cbind(simgrid, X)
  ## X may contain negative values, so it needs to be re-scaled
  weights <- (X+abs(min(X)))/max(X+abs(min(X)))
  indices <- sample(length(weighted.pts$X),
                    n.patch,
                    replace=FALSE,
                    prob=weights) 
  cbind(x=weighted.pts[indices,1],
        y=weighted.pts[indices,2])
}     

## landscape functions
locations.clumped.v2 <- function(size, acl, n.patch) {

  simgrid <- expand.grid(x=seq(1, size, by=1),
                         y=seq(1, size, by=1))

  landscape <- cn_2D(acl=acl, n=size, amp=1)$cn
  landscape <- landscape - min(landscape)
  landscape <- landscape/max(landscape)

  ## convert to landscape into habitable/non-habitable (binary) for
  ## our purposes
 # if(CTS) return(landscape)
  
  habitable <- landscape*0
  if(n.patch==0) return(habitable)
  habitable[order(landscape, decreasing=TRUE)[1:n.patch]] <- 1
  habitable
}     

## LKM: created this function (using your code) which generates the
## landscape using the above locations.clumped function
make.landscape <- function(ltc.size,
                           acl=acl,
                           frac.suitable) {

  ## ------------------------------------------------------------
  ## number of patches of each type
  n.patch <- frac.suitable*rep(ltc.size*ltc.size, 2)
  
  h1 <- locations.clumped.v2(ltc.size, acl[1], n.patch[1]) 
  h2 <- locations.clumped.v2(ltc.size, acl[2], n.patch[2])

  abind(h1, h2, along=3)
}
