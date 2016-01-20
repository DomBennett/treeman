# TODO: sort documentation

# @name get_Name
getNodeName <- function(tree, name, id) {
  tree@nodelist[[id]][[name]]
}

getNodesName <- function(tree, name, ids) {
  .get <- function(i) {
    getNodeName(tree, name, ids[i])
  }
  sapply(1:length(ids), .get)
}

# @name get_Children

getNodeChildren <- function(tree, id) {
  node <- tree@nodelist[[id]]
  node$children
}

getNodesChildren <- function(tree, ids) {
  sapply(ids, getNodeChildren, tree=tree)
}

# @name get_Age
#TODO: how to effectively handle unrooted trees, age has no meaning
getNodeAge <- function(tree, id) {
  node <- tree@nodelist[[id]]
  tree@age - node$predist
}

getNodesAge <- function(tree, ids) {
  ages <- sapply(ids, getNodeAge, tree=tree)
  data.frame(node=ids, age=ages, row.names=NULL)
}

getEdgeAge <- function(tree, id) {
  max <- getNodeAge(tree, tree@nodelist[[id]]$pre)
  min <- getNodeAge(tree, id)
  data.frame(edge=id, max, min)
}

getEdgesAge <- function(tree, ids) {
  maxs <- sapply(ids, function(tree, id) {
    getNodeAge(tree, tree@nodelist[[id]]$pre)
  }, tree=tree)
  mins <- sapply(ids, getNodeAge, tree=tree)
  data.frame(edge=ids, max=maxs, min=mins, row.names=NULL)
}

# @name getParent
getParent <- function(tree, ids) {
  prids <- getNodesPre(tree, ids)
  rf <- prids[[1]]
  mn_rnk <- 0
  for(n in prids[-1]) {
    rnk <- min(match(n, rf), na.rm=TRUE)
    if(rnk > mn_rnk) mn_rnk <- rnk
  }
  rf[mn_rnk]
}

# @name getPath
getPath <- function(tree, from, to) {
  pre_1 <- getNodePre(tree, from)
  pre_2 <- getNodePre(tree, to)
  parent <- pre_1[which(pre_1 %in% pre_2)[1]]
  path_1 <- c(from ,pre_1[!pre_1 %in% pre_2])
  path_2 <- c(pre_2[!pre_2 %in% pre_1], to)
  c(path_1, parent, path_2)
}

# @name get_Pre
getNodePre <- function(tree, id) {
  .get <- function(nd, prids) {
    prid <- tree@nodelist[[nd]]$pre
    if(!is.null(prid)) {
      prids <- c(prid, .get(prid, prids))
    }
    prids
  }
  .get(id, NULL)
}

getNodesPre <- function(tree, ids) {
  sapply(ids, getNodePre, tree=tree, simplify=FALSE)
}

# @name get_Lineage
getNodeLineage <- function(tree, id) {
  prids <- getNodePre(tree, id)
  lineage <- sapply(prids, function(n) tree@nodelist[[n]]$taxonym)
  if(length(lineage) > 0) {
    lineage <- c(tree@nodelist[[id]]$taxonym, lineage)
    lineage <- lineage[length(lineage):1]
    lineage <- unique(lineage)
  } else {
    lineage <- NULL
  }
  lineage
}

getNodesLineage <- function(tree, ids) {
  sapply(ids, getNodeLineage, tree=tree)
}

# @name get_Post
getNodePost <- function(tree, id) {
  .get <- function(nds, pstids) {
    new_nds <- c()
    for(nd in nds) {
      new_nds <- c(new_nds, tree@nodelist[[nd]]$post)
    }
    pstids <- c(pstids, new_nds)
    if(length(new_nds) > 0) {
      pstids <- .get(nds=new_nds, pstids=pstids)
    }
    pstids
  }
  .get(nds=node, pstids=NULL)
}

getNodesPost <- function(tree, ids) {
  sapply(ids, getNodePost, tree=tree)
}

# @name getSubtree
getSubtree <- function(tree, id) {
  pstids <- getNodePost(tree, id)
  ndlst <- tree@nodelist[c(id, pstids)]
  nd_prdst <- ndlst[[id]]$predist
  ndlst <- lapply(ndlst, function(x) {
    x$predist <- x$predist - nd_prdst
    x
  })
  ndlst[[id]]$pre <- NULL
  ndlst[[id]]$span <- 0
  new_tree <- new('TreeMan', nodelist=ndlst, root=id)
  .update(new_tree)
}
