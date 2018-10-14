## load source functions
library('abind')
library('parallel')

source('src/patchy.R')
source('src/prms.R')
source('src/space.R')
source('src/step.dist.R')

source('src/create.pop.R')
source('src/disperse.R')
source('src/disturb.R')
source('src/mutate.R')
#source('src/alt.mutate.R')
source('src/growth.R')
source('src/recombination.R')
source('src/select.R')
source('src/loss.R')

source('src/next.gen.recol.R')
## source('src/next.gen.no.recol.R')
source('src/meta.functions.R')

source('src/plotting.R')

## load and return loaded object
load.local <- function(file) {
  v <- load(file)
  stopifnot(length(v==1))
  get(v)
}
