# TODO: calc imbalance, calc tree dists

calcPhyDv <- function(tree, tips) {
  prnds <- unlist(getNodesPrenodes(tree, tips))
  counts <- table(prnds)
  prnds <- names(counts)[counts < length(tips)]
  spans <- unlist(lapply(tree@nodelist[c(tips, prnds)],
                         function(x) x$span))
  sum(spans)
}

calcFrPrp <- function(tree, tips) {
  .share <- function(node) {
    span <- tree@nodelist[[node]]$span
    children <- getNodeChildren(tree, node)
    if(!is.null(children)) {
      n <- length(children)
    } else {
      n <- 1
    }
    span/n
  }
  .calc <- function(tip) {
    nodes <- c(tip, getNodePrenodes(tree, tip))
    shares <- unlist(sapply(nodes, .share))
    sum(shares)
  }
  sapply(tips, .calc)
}

calcDstMtrx <- function(tree) {
  nodes <- names(tree@nodelist)
  cmbs <- expand.grid(nodes, nodes, stringsAsFactors=FALSE)
  .getDist <- function(cmb) {
    if(cmb[1] == cmb[2]) {
      return(0)
    }
    path <- getPath(tree, from=cmb[1], to=cmb[2])
    path_spans <- unlist(lapply(tree@nodelist[path], function(x) x$span))
    sum(path_spans)
  }
  res <- apply(cmbs, 1, .getDist)
  res <- matrix(res, ncol=length(nodes))
  colnames(res) <- rownames(res) <- nodes
  res
}