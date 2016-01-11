# TODO: sort documentation
# @name get_Children

getNodeChildren <- function(tree, node) {
  node <- tree@nodelist[[node]]
  node$children
}

getNodesChildren <- function(tree, nodes='all') {
  if(nodes[1] == 'all') {
    nodes <- tree@nodes
  }
  sapply(nodes, getNodeChildren, tree=tree)
}

# @name get_Age
#TODO: how to effectively handle unrooted trees, age has no meaning
getNodeAge <- function(tree, node) {
  node <- tree@nodelist[[node]]
  tree@age - node$predist
}

getNodesAge <- function(tree, nodes='all') {
  if(nodes[1] == 'all') {
    nodes <- c(tree@nodes, tree@tips)
  }
  ages <- sapply(nodes, getNodeAge, tree=tree)
  data.frame(node=nodes, age=ages, row.names=NULL)
}

getEdgeAge <- function(tree, edge) {
  max <- getNodeAge(tree, tree@nodelist[[edge]]$prenode)
  min <- getNodeAge(tree, edge)
  data.frame(edge, max, min)
}

getEdgesAge <- function(tree, edges='all') {
  if(edges[1] == 'all') {
    edges <- c(tree@nodes, tree@tips)
    edges <- edges[which(edges != tree@root)]
  }
  maxs <- sapply(edges, function(tree, edge) {
    getNodeAge(tree, tree@nodelist[[edge]]$prenode)
  }, tree=tree)
  mins <- sapply(edges, getNodeAge, tree=tree)
  data.frame(edge=edges, max=maxs, min=mins, row.names=NULL)
}

# @name getParent
getParent <- function(tree, nodes) {
  prnds <- getNodesPrenodes(tree, nodes)
  rf <- prnds[[1]]
  mn_rnk <- 0
  for(n in prnds[-1]) {
    rnk <- min(match(n, rf), na.rm=TRUE)
    if(rnk > mn_rnk) mn_rnk <- rnk
  }
  rf[mn_rnk]
}

# @name getPath
getPath <- function(tree, node_1, node_2) {
  # find parent between two nodes
  # get all prenodes to parent
  # return vector of prenodes from node_1 to node_2
  prenodes_1 <- getNodePrenodes(tree, node_1)
  prenodes_2 <- getNodePrenodes(tree, node_2)
  parent <- prenodes_1[which(prenodes_1 %in% prenodes_2)[1]]
  path_1 <- c(node_1 ,prenodes_1[!prenodes_1 %in% prenodes_2])
  path_2 <- c(prenodes_2[!prenodes_2 %in% prenodes_1], node_2)
  c(path_1, parent, path_2)
}

# @name get_Prenodes
getNodePrenodes <- function(tree, node) {
  .get <- function(nd, prnds) {
    prnd <- tree@nodelist[[nd]]$prenode
    if(!is.null(prnd)) {
      prnds <- c(prnd, .get(prnd, prnds))
    }
    prnds
  }
  .get(node, NULL)
}

getNodesPrenodes <- function(tree, nodes='all') {
  if(nodes[1] == 'all') {
    nodes <- c(tree@nodes, tree@tips)
    nodes <- nodes[which(nodes != tree@root)]
  }
  sapply(nodes, getNodePrenodes, tree=tree)
}

# @name get_Lineage
getNodeLineage <- function(tree, node) {
  prnds <- getNodePrenodes(tree, node)
  lineage <- c(tree@nodelist[[node]]$taxonym,
               sapply(prnds, function(n) tree@nodelist[[n]]$taxonym))
  lineage <- lineage[length(lineage):1]
  unique(lineage)
}

getNodesLineage <- function(tree, nodes) {
  sapply(nodes, getNodeLineage, tree=tree)
}

getNodesLineage(tree, tips(tree))

# @name get_Postnodes
getNodePostnodes <- function(tree, node) {
  .get <- function(nds, pstnds) {
    new_nds <- c()
    for(nd in nds) {
      new_nds <- c(new_nds, tree@nodelist[[nd]]$postnode)
    }
    pstnds <- c(pstnds, new_nds)
    if(length(new_nds) > 0) {
      pstnds <- .get(nds=new_nds, pstnds=pstnds)
    }
    pstnds
  }
  .get(nds=node, pstnds=NULL)
}

getNodesPostnodes <- function(tree, nodes='all') {
  if(nodes[1] == 'all') {
    nodes <- c(tree@nodes, tree@tips)
    nodes <- nodes[which(nodes != tree@root)]
  }
  sapply(nodes, getNodePostnodes, tree=tree)
}

# @name getSubtree
getSubtree <- function(tree, node) {
  pstnds <- getNodePostnodes(tree, node)
  ndlst <- tree@nodelist[c(node, pstnds)]
  nd_prdst <- ndlst[[node]]$predist
  ndlst <- lapply(ndlst, function(x) {
    x$predist <- x$predist - nd_prdst
    x
  })
  ndlst[[node]]$prenode <- NULL
  ndlst[[node]]$span <- 0
  new_tree <- new('TreeMan', nodelist=ndlst, root=node)
  treeman:::.update(new_tree)
}
