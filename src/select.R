## selection
select <- function(prms, pop)
  pop1 <- array(rbinom(length(pop), pop, prms$real.hab.fit),
                dim=dim(pop))


