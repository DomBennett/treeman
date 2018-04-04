# prevent sudden crashes by stopping incorrect arguments going to C code
.preCtest <- function(...) {
  argg <- c(as.list(environment()), list(...))
  bool <- FALSE
  for(arg in argg) {
    if(!is.numeric(arg) | length(arg) < 1) {
      bool <- TRUE
    }
  }
  if(bool) {
    stop('1 or more arguments are of the wrong type')
  }
  NULL
}

# MULTIPLE NDS

.getNdsPtidsFrmLst <- function(ndlst, ids, prinds, parallel, progress) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  out <- plyr::mlply(.data=l_data, .fun=.getNdPtidsFrmLst, ndlst=ndlst,
                     prinds=prinds, .parallel=parallel, .progress=progress)
  names(out) <- attr(out, 'split_labels')[,1]
  res <- out[1:length(out)]
  res
}

.getNdsPridsFrmLst <- function(ndlst, ids, prinds,
                               parallel, progress) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  out <- plyr::mlply(.data=l_data, .fun=.getNdPridsFrmLst, ndlst=ndlst,
                     prinds=prinds, .parallel=parallel, .progress=progress)
  names(out) <- attr(out, 'split_labels')[,1]
  res <- out[1:length(out)]
  res
}

.getNdsPDFrmLst <- function(ndlst, ids, prinds,
                            parallel, progress) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  out <- plyr::mdply(.data=l_data, .fun=.getNdPDFrmLst, prinds=prinds,
                     ndlst=ndlst, .parallel=parallel, .progress=progress)
  res <- out[ ,2]
  names(res) <- out[ ,1]
  res
}

.getNdsKidsFrmLst <- function(ndlst, ids, prinds, tinds,
                              parallel, progress) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- plyr::mlply(.data=l_data, .fun=.getNdKidsFrmLst,
                     ndlst=ndlst, prinds=prinds,
                     tinds=tinds, .parallel=parallel, .progress=progress)
  names(res) <- ids
  res[1:length(res)]
}

.getNdsPrdstsFrmLst <- function(ndlst, ids, prinds,
                               parallel, progress) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  out <- plyr::mdply(.data=l_data, .fun=.getNdPrdstsFrmLst, prinds=prinds,
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
.getNdPridsFrmLst <- function(ndlst, prinds, id) {
  prid <- ndlst[[id]][['prid']]
  nids <- names(ndlst)
  prind <- which(nids == prid)
  .preCtest(prind, prinds)
  res <- .Call("cGetNdPrids", PACKAGE="treeman",
               as.integer(prind),
               as.integer(prinds))
  nids[res]
}

.getNdPrdstsFrmLst <- function(ndlst, prinds, id) {
  prids <- .getNdPridsFrmLst(ndlst, prinds, id)
  sum(vapply(ndlst[prids], function(x) x[['spn']], numeric(1))) +
    ndlst[[id]][['spn']]
}

#' @useDynLib treeman cGetNdPtids
.getNdPtidsFrmLst <- function(ndlst, prinds, id) {
  nids <- names(ndlst)
  id <- which(nids == id)
  .preCtest(id, prinds)
  res <- .Call("cGetNdPtids", PACKAGE="treeman",
               as.integer(id),
               as.integer(prinds))
  nids[which(res > 0)]
}

.getNdKidsFrmLst <- function(ndlst, prinds, tinds, id) {
  ptids <- .getNdPtidsFrmLst(ndlst, prinds, id)
  tids <- names(ndlst)[tinds]
  ptids[ptids %in% tids]
}

.getNdPDFrmLst <- function(ndlst, prinds, id) {
  ptids <- .getNdPtidsFrmLst(ndlst, prinds, id)
  if(length(ptids) > 0) {
    res <- sum(vapply(ndlst[ptids], function(x) x[['spn']], numeric(1)))
  } else {
    res <- 0
  }
  res
}

# TREE FUNCTIONS

.getTreeAgeFrmLst <- function(ndlst, prinds, tids, parallel) {
  prdsts <- .getNdsPrdstsFrmLst(ndlst, ids=tids, prinds=prinds,
                                parallel=parallel, progress="none")
  max(prdsts)
}

# SLOT

.getSltSpns <- function(ndlst) {
  .get <- function(x) {
    x[['spn']]
  }
  res <- vapply(ndlst, .get, numeric(1))
  names(res) <- NULL
  res
}

# INDS

.getPrinds <- function(ndlst) {
  # pre-node index
  .get <- function(x) {
    x[['prid']]
  }
  res <- vapply(ndlst, .get, character(1))
  match(res, names(ndlst))
}

.getTinds <- function(ndlst) {
  # tip-node index
  .get <- function(x) {
    length(x[['ptid']]) == 0
  }
  res <- vapply(ndlst, .get, logical(1))
  names(res) <- NULL
  which(res)
}

# SPECIAL

#' @useDynLib treeman cGetNdmtrx
.getNdmtrxFrmLst <- function(ndlst, shared=FALSE, ...) {
  # return matrix of 01s for ids that descend from 
  message('Note, trees with `ndmtrx` cannot be saved and loaded using `save()` or `savehistory()`.',
          ' Loading from these files may cause unusual behaviour.')
  prids <- vapply(ndlst, function(x) x[['prid']], character(1))
  nids <- names(prids)
  prids <- match(prids, nids)
  qry_ids <- 1:length(nids)
  .preCtest(length(nids), qry_ids, prids)
  res <- .Call("cGetNdmtrx", PACKAGE="treeman",
               as.integer(length(nids)),
               as.integer(qry_ids),
               as.integer(prids))
  res <- bigmemory::as.big.matrix(res, shared=shared, ...)
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