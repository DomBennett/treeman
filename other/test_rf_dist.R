# Compare RF dist between ape and treeman
# TODO: test with different numbers of tips

# Libs
library(treeman)

# Functions
mkphylo <- function(x) {
  ape::compute.brlen(ape::rtree(10))
}

# For
phylos <- lapply(1:10, mkphylo)
rf_dists_tm <- rf_dists_ape <- matrix(0, nrow = length(phylos),
                                      ncol = length(phylos))
for (i in seq_along(phylos)) {
  cat(i, "\n", sep = "")
  for (j in seq_along(phylos)) {
    #dst <- suppressWarnings(ape::dist.topo(x = phylos[[i]], y = phylos[[j]]))
    dst <- ape::dist.topo(x = ape::unroot(phylos[[i]]),
                          y = ape::unroot(phylos[[j]]))
    rf_dists_ape[i, j] <- dst
    dst <- calcDstRF(as(phylos[[i]], 'TreeMan'),
                     as(phylos[[j]], 'TreeMan'))
    rf_dists_tm[i, j] <- dst
  }
}

which(rf_dists_ape != rf_dists_tm)
rf_dists_ape == rf_dists_tm
par(mfrow=c(1, 2))

i <- 1
j <- 2
plot(phylos[[i]])
plot(phylos[[j]])
rf_dists_ape[i, j]
rf_dists_tm[i, j]
tree_1 <- as(phylos[[i]], 'TreeMan')
tree_2 <- as(phylos[[j]], 'TreeMan')


