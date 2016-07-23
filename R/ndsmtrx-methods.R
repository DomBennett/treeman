
# MULTIPLE NDS
.getNdsPrdstsFrmMtrx <- function(ndmtrx, all_ids, ids, spns,
                                 parallel, progress) {
  .get <- function(x) {
    sum(spns[x])
  }
  res <- plyr::adply(ndmtrx[ ,all_ids %in% ids], .margins=2,
                     .fun=.get, .parallel=parallel, .progress=progress)[ ,2]
  res <- res + spns[all_ids %in% ids]
  names(res) <- ids
  res
}

.getNdsKidsFrmMtrx <- function(ndmtrx, all_ids, ids, tids,
                               parallel, progress) {
  .get <- function(x) {
    all_ids[x & all_ids %in% tids]
  }
  res <- plyr::alply(ndmtrx[all_ids %in% ids, ], .margins=1,
                     .fun=.get, .parallel=parallel, .progress=progress)
  names(res) <- ids
  res <- res[1:length(res)]
  res
}

.getNdsPtidsFrmMtrx <- function(ndmtrx, all_ids, ids,
                                parallel, progress) {
  .get <- function(x) {
    all_ids[x]
  }
  res <- plyr::alply(ndmtrx[all_ids %in% ids, ], .margins=1,
                     .fun=.get, .parallel=parallel, .progress=progress)
  names(res) <- ids
  res <- res[1:length(res)]
  res
}

.getNdsPridsFrmMtrx <- function(ndmtrx, all_ids, ids,
                                parallel, progress) {
  .get <- function(x) {
    all_ids[x]
  }
  res <- plyr::alply(ndmtrx[ ,all_ids %in% ids], .margins=2,
                     .fun=.get, .parallel=parallel, .progress=progress)
  names(res) <- ids
  res <- res[1:length(res)]
  res
}

.getNdsPDFrmMtrx <- function(ndmtrx, all_ids, ids, spns,
                             parallel, progress) {
  .get <- function(x) {
    sum(spns[x])
  }
  res <- plyr::adply(ndmtrx[all_ids %in% ids, ], .margins=1,
                     .fun=.get, .parallel=parallel, .progress=progress)[ ,2]
  names(res) <- ids
  res <- res[1:length(res)]
  res
}

# TREE

.getTreeAgeFrmMtrx <- function(ndmtrx, all_ids, tids, spns, parallel) {
  res <- .getNdsPrdstsFrmMtrx(ndmtrx, all_ids, tids, spns,
                              parallel=parallel, progress='none')
  max(res)
}
