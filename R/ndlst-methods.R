
# MULTIPLE NDS

.getNdsPtidsFrmLst <- function(ndlst, ids, parallel, progress) {
  prids <- .getSltPrids(ndlst, parallel=parallel)
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  out <- plyr::mlply(.data=l_data, .fun=.getNdPtidFrmLst, ndlst=ndlst,
                     prids=prids, .parallel=parallel, .progress=progress)
  names(out) <- attr(out, 'split_labels')[,1]
  res <- out[1:length(out)]
  res
}

.getNdsPridsFrmLst <- function(ndlst, ids,
                               parallel, progress) {
  prids <- .getSltPrids(ndlst, parallel=parallel)
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  out <- plyr::mlply(.data=l_data, .fun=.getNdPridsFrmLst, ndlst=ndlst,
                     prids=prids, .parallel=parallel, .progress=progress)
  names(out) <- attr(out, 'split_labels')[,1]
  res <- out[1:length(out)]
  res
}

.getNdsPDFrmLst <- function(ndlst, ids,
                            parallel, progress) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  out <- plyr::mdply(.data=l_data, .fun=.getNdPDFrmLst,
                     ndlst=ndlst, .parallel=parallel, .progress=progress)
  res <- out[ ,2]
  names(res) <- out[ ,1]
  res
}

.getNdsKidsFrmLst <- function(ndlst, ids,
                              parallel, progress) {
  tids <- .getSltTids(ndlst, parallel)
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- plyr::mlply(.data=l_data, .fun=.getNdKidsFrmLst, ndlst=ndlst,
                     tids=tids, .parallel=parallel, .progress=progress)
  names(res) <- ids
  res[1:length(res)]
}

.getNdsPrdstsFrmLst <- function(ndlst, ids,
                               parallel, progress) {
  prids <- .getSltPrids(ndlst, parallel)
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  out <- plyr::mdply(.data=l_data, .fun=.getNdPrdstsFrmLst, prids=prids,
                     ndlst=ndlst, .parallel=parallel, .progress=progress)
  res <- out[ ,2]
  names(res) <- out[ ,1]
  res
}

# SINGLE ND

.getNdSstrFrmLst <- function(ndlst, id) {
  prid <- ndlst[[id]][['prid']][[1]]
  ptids <- ndlst[[prid]][['ptid']]
  ptids[ptids != id]
}

#' @useDynLib treeman cGetNdPrids
.getNdPridsFrmLst <- function(ndlst, prids, id) {
  prid <- ndlst[[id]][['prid']]
  nids <- names(ndlst)
  prid <- which(nids == prid)
  prids <- match(prids, nids)
  res <- .Call("cGetNdPrids", PACKAGE="treeman",
               as.integer(prid),
               as.integer(prids))
  nids[res]
}

.getNdPrdstsFrmLst <- function(ndlst, prids, id) {
  prids <- .getNdPridsFrmLst(ndlst, prids, id)
  sum(sapply(ndlst[prids], function(x) x[['spn']])) +
    ndlst[[id]][['spn']]
}

#' @useDynLib treeman cGetNdPtids
.getNdPtidsFrmLst <- function(ndlst, prids, id) {
  nids <- names(ndlst)
  id <- which(nids == id)
  prids <- match(prids, nids)
  res <- .Call("cGetNdPtids", PACKAGE="treeman",
               as.integer(id),
               as.integer(prids))
  nids[which(res > 0)]
}

.getNdKidsFrmLst <- function(ndlst, prids, tids, id) {
  ptids <- .getNdPtidsFrmLst(ndlst, prids, id)
  ptids[ptids %in% tids]
}

.getNdPDFrmLst <- function(ndlst, prids, id) {
  ptids <- .getNdPtidsFrmLst(ndlst, prids, id)
  if(length(ptids) > 0) {
    res <- sum(sapply(ndlst[ptids], function(x) x[['spn']]))
  } else {
    res <- 0
  }
  res
}

# TREE FUNCTIONS

.getTreeAgeFrmLst <- function(ndlst, parallel) {
  tids <- .getSltTids(ndlst, parallel)
  prdsts <- .getNdsPrdstsFrmLst(ndlst, tids,
                               parallel=parallel, progress="none")
  max(prdsts)
}

# SLOT

.getSltPrids <- function(ndlst, parallel) {
  .get <- function(x) {
    x[['prid']]
  }
  res <- plyr::ldply(ndlst, .fun=.get, .parallel=parallel)
  res[ ,2]
}

.getSltSpns <- function(ndlst, parallel) {
  .get <- function(x) {
    x[['spn']]
  }
  res <- plyr::ldply(ndlst, .fun=.get, .parallel=parallel)
  res[ ,2]
}

.getSltTids <- function(ndlst, parallel) {
  .get <- function(x) {
    length(x[['ptid']]) == 0
  }
  res <- plyr::ldply(ndlst, .fun=.get, .parallel=parallel)
  res[ ,1][res[ ,2]]
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