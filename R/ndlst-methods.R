# Attemp for making getNdsMat run in parallel
# ... actually made it slower
# ntids <- length(tids)
# n <- foreach::getDoParWorkers()
# nparts <- ntids %/% n
# parts <- c(seq(1, ntids - 1, nparts), ntids + 1)
# res <- foreach (i=2:length(parts), .combine="cbind") %dopar% {
#   tids <- tids[parts[i-1]:(parts[i] - 1)]
#   res <- .Call("cGetNdsMat", PACKAGE="treeman",
#                as.integer(length(nids)),
#                as.integer(tids),
#                as.integer(prids))
#   res
# }

# MULTIPLE NDS
#' @useDynLib treeman cGetNdsMat
.getNdsMat <- function(ndlst, qry_ids) {
  # return matrix of 01s for ids that descend from 
  prids <- sapply(ndlst, function(x) x[['prid']])
  nids <- names(prids)
  prids <- match(prids, nids)
  tids <- match(qry_ids, nids)
  res <- .Call("cGetNdsMat", PACKAGE="treeman",
               as.integer(length(nids)),
               as.integer(tids),
               as.integer(prids))
  res <- res > 0
  res
}

# SINGLE ND

.getNdSstr <- function(ndlst, id) {
  prid <- ndlst[[id]][['prid']][[1]]
  ptids <- ndlst[[prid]][['ptid']]
  ptids[ptids != id]
}

#' @useDynLib treeman cGetNdPrids
.getNdPrids <- function(ndlst, id) {
  prids <- sapply(ndlst, function(x) x[['prid']])
  prid <- ndlst[[id]][['prid']]
  nids <- names(prids)
  prid <- which(nids == prid)
  prids <- match(prids, nids)
  res <- .Call("cGetNdPrids", PACKAGE="treeman",
               as.integer(prid),
               as.integer(prids))
  res
}

.getNdPrdst <- function(ndlst, id) {
  prids <- .getNdPrids(ndlst, id)
  sum(sapply(ndlst[prids], function(x) x[['spn']])) +
    ndlst[[id]][['spn']]
}

#' @useDynLib treeman cGetNdPtids
.getNdPtids <- function(ndlst, id) {
  prids <- sapply(ndlst, function(x) x[['prid']])
  nids <- names(prids)
  id <- which(nids == id)
  prids <- match(prids, nids)
  res <- .Call("cGetNdPtids", PACKAGE="treeman",
               as.integer(id),
               as.integer(prids))
  which(res > 0)
}

.getNdKids <- function(ndlst, id) {
  ptids <- .getNdPtids(ndlst, id)
  kids <- sapply(ndlst[ptids], function(x) length(x[['ptid']]) == 0)
  ptids[as.logical(kids)]
}

.getNdPD <- function(ndlst, id) {
  ptids <- .getNdPtids(ndlst, id)
  if(length(ptids) > 0) {
    res <- sum(sapply(ndlst[ptids], function(x) x[['spn']]))
  } else {
    res <- 0
  }
  res
}

# TREE FUNCTIONS

.getTreeAge <- function(ndlst) {
  tids <- sapply(ndlst, function(x) length(x[['ptid']]) == 0)
  tids <- as.integer(which(tids))
  tip_prdsts <- sapply(tids, .getNdPrdst, ndlst=ndlst)
  max(tip_prdsts)
}