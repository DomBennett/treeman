# TODO: calc evolutionary distinctness, calc imbalance, calc tree dists

calcPhyDiv <- function(tree, nodes) {
  #TODO
  # 1. Get parent
  # 2. Get postnodes
  # 3. sum the branch lengths of the postnodes
}

calcFairProp <- function(tree, tips='all') {
  .share <- function(node) {
    span <- tree@nodelist[[node]]$span
    nchildren <- length(getNodeChildren(tree, node))
    span/nchildren
  }
  .calc <- function(node) {
    prnds <- getNodePrenodes(tree, node)
    shares <- unlist(sapply(prnds, .share))
    sum(shares)
  }
  if(tips == 'all') {
    tips <- tree@tips
  }
  sapply(tips, .calc)
}

calcDstMtrx <- function(tree) {
  cmbs <- expand.grid(nodes, nodes, stringsAsFactors=FALSE)
  .getDist <- function(cmb) {
    if(cmb[1] == cmb[2]) {
      return(0)
    }
    path <- getPath(tree, node_1=cmb[1], node_2=cmb[2])
    path_spans <- unlist(lapply(tree@nodelist[path], function(x) x$span))
    sum(path_spans)
  }
  res <- apply(cmbs, 1, .getDist)
  res <- matrix(res, ncol=length(nodes))
  colnames(res) <- rownames(res) <- nodes
  res
}