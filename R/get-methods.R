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
getPath <- function(tree, from, to) {
  prenodes_1 <- getNodePrenodes(tree, from)
  prenodes_2 <- getNodePrenodes(tree, to)
  parent <- prenodes_1[which(prenodes_1 %in% prenodes_2)[1]]
  path_1 <- c(from ,prenodes_1[!prenodes_1 %in% prenodes_2])
  path_2 <- c(prenodes_2[!prenodes_2 %in% prenodes_1], to)
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

getNodesPrenodes <- function(tree, nodes) {
  sapply(nodes, getNodePrenodes, tree=tree, simplify=FALSE)
}

# @name get_Lineage
getNodeLineage <- function(tree, node) {
  prnds <- getNodePrenodes(tree, node)
  lineage <- sapply(prnds, function(n) tree@nodelist[[n]]$taxonym)
  if(length(lineage) > 0) {
    lineage <- c(tree@nodelist[[node]]$taxonym, lineage)
    lineage <- lineage[length(lineage):1]
    lineage <- unique(lineage)
  } else {
    lineage <- NULL
  }
  lineage
}

getNodesLineage <- function(tree, nodes) {
  sapply(nodes, getNodeLineage, tree=tree)
}

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

getNodesPostnodes <- function(tree, nodes) {
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
  .update(new_tree)
}
