## Function to generate autocorrelated noise in 2D

## Function producing white noise with expected variance sigma^2

## Gaussian white noise
wn_gauss <- function(n=1,sigma=1)
  return(rnorm(n,mean=0,sd=sigma))

## locations in (0,1) of the n points
xvec <- function(n,dx=NULL) {
  if (is.null(dx))
    return((1:n-0.5)/n)
  else
    return((1:n-0.5)*dx)
}

## Filter function
##
## l0 is a "magic" multiple of the autocorrelation length where the
## function value has dropped to mostly zero; to preserve the white
## noise amplitude, the filter is normalized so that sum(f[i]^2)==1
## (the kern() functions as given are normalized to one, which is of
## course unnecessary given the eventual renormalization of f) 2D
## integration with adaptIntegrate() from package "cubature"
require("cubature")

## Gaussian filter
flt_gauss_2D <- function(dx, acl=1) {
  l0 <- 3.5
  m <- ceiling(l0*acl/dx)
  x <- (1:m)*dx
  kern <- function(y) exp(-sum(y^2)/(2*acl^2))/(2*pi*acl^2)
  r <- matrix(nrow=m,ncol=m,byrow=FALSE)
  for(j in 1:m)
    for(i in 1:m)
      r[i,j] =
        adaptIntegrate(kern,c(x[i]-dx,x[j]-dx),c(x[i],x[j]))$integral
  f <- rbind(cbind(r[m:1,m:1],r[m:1,]), cbind(r[,m:1],r))
  return(f/sqrt(sum(f^2)))
}

## generate correlated noise
## n     ... number of points equidistant in (0,1)
## acl   ... desired autocorrelation length
## amp   ... amplitude of noise
## wn    ... function generating noise vector taking args (n,amp)
## flt   ... function generating impulse response vector taking args (x,acl)
##                               where x is a vector of locations

## for 2D
cn_2D <- function(acl,n=16,amp=1,wn=wn_gauss,flt=flt_gauss_2D) {
  print(paste("Will generate a",n,"by",n,"grid."))
  f <- flt(1/n,acl) # m x m-matrix of filter response
  m <- dim(f)[1]
  print(paste("Generated a",m,"by",m,"filter."))
  if (m>n) {
    print("WARNING: Filter too long.")
    ## truncate filter
    f <- f[(1:n)+(m-n)/2,]
  }
  ## pad filter with zeros so that dim==(n,nx)
  if (m<n) {
    f <- cbind(matrix(0,m,(n-m)/2),f,matrix(0,m,(n-m)/2))
    f <- rbind(matrix(0,(n-m)/2,n),f,matrix(0,(n-m)/2,n))
  }
  else {
    f <- cbind(matrix(0,n,(n-m)/2),f,matrix(0,n,(n-m)/2))
  }
  ## now create n x n uncorrelated noise landscape
  w <- matrix(wn(n*n,amp),n,n) # n x nx matrix of noise peaks
  h <- Re(fft(fft(f)*fft(w),inverse=TRUE))/(n*n)
  return(list(cn=h,flt=f,n=n,m=m))
}

## bilinear interpolation
## returns value at z by inter- and extrapolating on f
## length(z)==2, dim(f)=c(n,n)
bilerp <- function(z,f) {
  n <- sqrt(length(f))
  x <- xvec(n) # all points in the interior of (0,1), x[i+1]-x[i]=1/n
  j <- max(min(ceiling((z[1]-x[1])*n),n-1),1) # find index, assuming dx == 1/n
  k <- max(min(ceiling((z[2]-x[1])*n),n-1),1) # keeping index in range 1:(n-1)
  y1 <- f[k,j]
  y2 <- f[k,j+1]
  y3 <- f[k+1,j+1]
  y4 <- f[k+1,j]
  
  t <- (z[1]-x[j])*n
  u <- (z[2]-x[k])*n
  
  return((1-t)*(1-u)*y1+t*(1-u)*y2+t*u*y3+(1-t)*u*y4)
}

resample_landscape <- function(u,m=ncol(u)) {
  if (m==ncol(u)) return(u)
  
  x <- xvec(m)
  r <- matrix(NA,m,m)
  for (j in 1:m) {
    for (i in 1:m) {
      r[j,i] <- bilerp(c(x[i],x[j]),u)
    }
  }
  return(r)
}
