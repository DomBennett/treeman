#' @name get_Children
#' @title 
#' @description 
#' @details 
#' @param
#' @seealso
#' @export
#' @examples

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

#' @name get_Age
#' @title 
#' @description 
#' @details 
#' @param
#' @seealso
#' @export
#' @examples
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
  age_1 <- getNodeAge(tree, tree@nodelist[[edge]]$prenode)
  age_2 <- getNodeAge(tree, edge)
  data.frame(edge, age_1, age_2)
}

getEdgesAge <- function(tree, edges='all') {
  if(edges[1] == 'all') {
    edges <- c(tree@nodes, tree@tips)
    edges <- edges[which(edges != tree@root)]
  }
  age_1s <- sapply(edges, function(tree, edge) {
    getNodeAge(tree, tree@nodelist[[edge]]$prenode)
  }, tree=tree)
  age_2s <- sapply(edges, getNodeAge, tree=tree)
  data.frame(edge=edges, age_1=age_1s, age_2=age_2s, row.names=NULL)
}

#' @name getParent
#' @title 
#' @description 
#' @details 
#' @param
#' @seealso
#' @export
#' @examples

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

#' @name get_Prenodes
#' @title 
#' @description 
#' @details 
#' @param
#' @seealso
#' @export
#' @examples

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

#' @name get_Postnodes
#' @title 
#' @description 
#' @details 
#' @param
#' @seealso
#' @export
#' @examples

getNodePostnodes <- function() {
  #TODO
}

getNodesPostnodes <- function() {
  #TODO
}