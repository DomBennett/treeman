
# MULTIPLE NDS

.getNdsKidsFrmLst <- function(ndlst, ids,
                              parallel=FALSE,
                              progress="none") {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- plyr::mlply(.data=l_data, .fun=.getNdKidsFrmLst,
                     ndlst=ndlst, .parallel=parallel, .progress=progress)
  names(res) <- ids
  res[1:length(res)]
}

# SINGLE ND

.getNdSstrFrmLst <- function(ndlst, id) {
  prid <- ndlst[[id]][['prid']][[1]]
  ptids <- ndlst[[prid]][['ptid']]
  ptids[ptids != id]
}

#' @useDynLib treeman cGetNdPrids
.getNdPridsFrmLst <- function(ndlst, id) {
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

.getNdPrdstFrmLst <- function(ndlst, id) {
  prids <- .getNdPridsFrmLst(ndlst, id)
  sum(sapply(ndlst[prids], function(x) x[['spn']])) +
    ndlst[[id]][['spn']]
}

#' @useDynLib treeman cGetNdPtids
.getNdPtidsFrmLst <- function(ndlst, id) {
  prids <- sapply(ndlst, function(x) x[['prid']])
  nids <- names(prids)
  id <- which(nids == id)
  prids <- match(prids, nids)
  res <- .Call("cGetNdPtids", PACKAGE="treeman",
               as.integer(id),
               as.integer(prids))
  which(res > 0)
}

.getNdKidsFrmLst <- function(ndlst, id) {
  ptids <- .getNdPtidsFrmLst(ndlst, id)
  kids <- sapply(ndlst[ptids], function(x) length(x[['ptid']]) == 0)
  names(ndlst)[ptids[as.logical(kids)]]
}

.getNdPDFrmLst <- function(ndlst, id) {
  ptids <- .getNdPtidsFrmLst(ndlst, id)
  if(length(ptids) > 0) {
    res <- sum(sapply(ndlst[ptids], function(x) x[['spn']]))
  } else {
    res <- 0
  }
  res
}

# TREE FUNCTIONS

.getTreeAgeFrmLst <- function(ndlst) {
  tids <- sapply(ndlst, function(x) length(x[['ptid']]) == 0)
  tids <- as.integer(which(tids))
  tip_prdsts <- sapply(tids, .getNdPrdstFrmLst, ndlst=ndlst)
  max(tip_prdsts)
}

# SPECIAL

#' @useDynLib treeman cGetNdmtrx
.getNdmtrxFrmLst <- function(ndlst) {
  # return matrix of 01s for ids that descend from 
  prids <- sapply(ndlst, function(x) x[['prid']])
  nids <- names(prids)
  prids <- match(prids, nids)
  qry_ids <- 1:length(nids)
  res <- .Call("cGetNdmtrx", PACKAGE="treeman",
               as.integer(length(nids)),
               as.integer(qry_ids),
               as.integer(prids))
  res <- res > 0
  res
}
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