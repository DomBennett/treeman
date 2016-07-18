#' @useDynLib treeman getPrids
.getPrids <- function(ndlst, id) {
  prids <- sapply(ndlst, function(x) x[['prid']])
  prid <- ndlst[[id]][['prid']]
  nids <- names(prids)
  prid <- which(nids == prid)
  prids <- match(prids, nids)
  res <- .Call("getPrids", PACKAGE="treeman",
               as.integer(prid),
               as.integer(prids))
  res
}

.getPrdst <- function(ndlst, id) {
  prids <- .getPrids(ndlst, id)
  sum(sapply(ndlst[prids], function(x) x[['spn']])) +
    ndlst[[id]][['spn']]
}

#' @useDynLib treeman getPtids
.getPtids <- function(ndlst, id) {
  prids <- sapply(ndlst, function(x) x[['prid']])
  nids <- names(prids)
  id <- which(nids == id)
  prids <- match(prids, nids)
  res <- .Call("getPtids", PACKAGE="treeman",
               as.integer(id),
               as.integer(prids))
  which(res > 0)
}

.getKids <- function(ndlst, id) {
  ptids <- .getPtids(ndlst, id)
  kids <- sapply(ndlst[ptids], function(x) length(x[['ptid']]) == 0)
  ptids[as.logical(kids)]
}

.getPD <- function(ndlst, id) {
  ptids <- .getPtids(ndlst, id)
  sum(sapply(ndlst[ptids], function(x) x[['spn']]))
}