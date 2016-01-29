# TODO: calc imbalance

calcDstTrp <- function(tree_1, tree_2, nrmlsd=FALSE) {
  .count <- function(i) {
    o1 <- getOutgroup(tree_1, cmbs[ ,i])
    o2 <- getOutgroup(tree_2, cmbs[ ,i])
    if (length(o1) != length(o2) || o1 != o2) {
      cntr <<- cntr + 1
    }
    NULL
  }
  shrd <- tree_1@tips[tree_1@tips %in% tree_2@tips]
  cmbs <- combn(shrd, 3)
  cntr <- 0
  sapply(1:ncol(cmbs), .count)
  if (nrmlsd) {
    cntr <- cntr/ncol(cmbs)
  }
  cntr
}

calcOvrlp <- function(tree, ids_1, ids_2, nrmlsd=FALSE) {
  spans <- getNodesSlot(tree, name='span', tree@all)
  names(spans) <- tree@all
  ids_1 <- c(unique(unlist(getNodesPrid(tree, ids_1))),
             ids_1)
  ids_2 <- c(unique(unlist(getNodesPrid(tree, ids_2))),
             ids_2)
  ovrlp <- sum(spans[ids_2[ids_2 %in% ids_1]])
  if(nrmlsd) {
    ovrlp <- ovrlp/tree@pd
  }
  ovrlp
}

calcDstBLD <- function(tree_1, tree_2, nrmlsd=FALSE) {
  n1 <- tree_1@nodes[!tree_1@nodes == tree_1@root]
  n2 <- tree_2@nodes[!tree_2@nodes == tree_2@root]
  c1 <- getNodesChildren(tree_1, n1)
  c2 <- getNodesChildren(tree_2, n2)
  s1 <- getNodesSlot(tree_1, name="span", ids=n1)
  s2 <- getNodesSlot(tree_2, name="span", ids=n2)
  d1 <- s2[match(c1, c2)]
  d1[which(is.na(d1))] <- 0
  d1 <- s1 - d1
  d2 <- s1[match(c2, c1)]
  d2[which(is.na(d2))] <- 0
  d2 <- s2 - d2
  d <- sqrt(sum(c(d1^2, d2^2)))
  if(nrmlsd) {
    max_d <- sqrt(sum(c(s1^2, s2^2)))
    d <- d/max_d
  }
  d
}

calcDstRF <- function(tree_1, tree_2, nrmlsd=FALSE) {
  n1 <- tree_1@nodes[!tree_1@nodes == tree_1@root]
  n2 <- tree_2@nodes[!tree_2@nodes == tree_2@root]
  c1 <- getNodesChildren(tree_1, n1)
  c2 <- getNodesChildren(tree_2, n2)
  d <- sum(!c1 %in% c2) + sum(!c2 %in% c1)
  if(nrmlsd) {
    max_d <- (length(n1) + length(n2))
    d <- d/max_d
  }
  d
}

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