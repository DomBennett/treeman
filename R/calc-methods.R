# TODO: calc imbalance, calc tree dists

calcPhyDv <- function(tree, ids) {
  prids <- unlist(getNodesPrid(tree, ids))
  counts <- table(prids)
  prids <- names(counts)[counts < length(ids)]
  spans <- unlist(lapply(tree@nodelist[c(ids, prids)],
                         function(n) n[['span']]))
  sum(spans)
}

calcFrPrp <- function(tree, ids) {
  .share <- function(id) {
    span <- tree@nodelist[[id]][['span']]
    children <- getNodeChildren(tree, id)
    if(!is.null(children)) {
      n <- length(children)
    } else {
      n <- 1
    }
    span/n
  }
  .calc <- function(tip) {
    ids <- c(tip, getNodePrid(tree, tip))
    shares <- unlist(sapply(ids, .share))
    sum(shares)
  }
  sapply(ids, .calc)
}

calcDstMtrx <- function(tree, ids) {
  cmbs <- expand.grid(ids, ids, stringsAsFactors=FALSE)
  .getDist <- function(cmb) {
    if(cmb[1] == cmb[2]) {
      return(0)
    }
    path <- getPath(tree, from=cmb[1], to=cmb[2])
    path_spans <- unlist(lapply(tree@nodelist[path], function(n) n[['span']]))
    sum(path_spans)
  }
  res <- apply(cmbs, 1, .getDist)
  res <- matrix(res, ncol=length(ids))
  colnames(res) <- rownames(res) <- ids
  res
}