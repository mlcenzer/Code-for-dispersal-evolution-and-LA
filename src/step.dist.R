make.step.dist<-function(max.dist, n.bins){
  ## create distance class object for possible dispersal steps
  step.grid <- expand.grid(x=-max.dist:max.dist,
                           y=-max.dist:max.dist)
  ## calculate euclidean step distances
  step.dist <- cbind(step.grid, dist=sqrt(rowSums(step.grid^2)))

  step.dist <- step.dist[step.dist[,'dist']<max.dist,]
  step.dist <- step.dist[step.dist$dist>0,] ## drop 0th bin
  
  ## add distance class column
  bin <- ceiling(step.dist[,'dist']/(max.dist/n.bins))
  step.dist <- cbind(step.dist, class=bin)

  dist.class <- split(1:nrow(step.dist), step.dist[,'class'])
  if(length(dist.class)<n.bins) {
    cat('ERROR: dist.class missing element(s)\n')
    break
  }
  list(step.dist, dist.class)
}
