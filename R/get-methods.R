# TODO: sort documentation

getOutgroup <- function(tree, ids) {
  .cntr <- function(id) {
    children <- tree@nodelist[[id]][['children']]
    sum(ids %in% children)
  }
  prnt <- getParent(tree, ids)
  ptids <- tree@nodelist[[prnt]][['ptid']]
  cnts <- sapply(ptids, .cntr)
  outnd <- names(cnts)[which.min(cnts)]
  children <- tree@nodelist[[outnd]][['children']]
  if(length(children) == 0) {
    return(outnd)
  }
  outgroup <- ids[ids %in% children]
  if(length(outgroup) > 1) {
    return(NULL)
  }
  outgroup
}

# @name get_Name
getNodeSlot <- function(tree, name, id) {
  tree@nodelist[[id]][[name]]
}

getNodesSlot <- function(tree, name, ids) {
  .get <- function(i) {
    getNodeSlot(tree, name, ids[i])
  }
  sapply(1:length(ids), .get)
}

# @name get_Children

getNodeChildren <- function(tree, id) {
  node <- tree@nodelist[[id]]
  node[['children']]
}

getNodesChildren <- function(tree, ids) {
  sapply(ids, getNodeChildren, tree=tree, simplify=FALSE)
}

# @name get_Age
#TODO: how to effectively handle unrooted trees, age has no meaning
getNodeAge <- function(tree, id) {
  node <- tree@nodelist[[id]]
  tree@age - node[['prdst']]
}

getNodesAge <- function(tree, ids) {
  ages <- sapply(ids, getNodeAge, tree=tree)
  data.frame(node=ids, age=ages, row.names=NULL)
}

getEdgeAge <- function(tree, id) {
  max <- getNodeAge(tree, tree@nodelist[[id]][['prid']])
  min <- getNodeAge(tree, id)
  data.frame(edge=id, max, min)
}

getEdgesAge <- function(tree, ids) {
  maxs <- sapply(ids, function(tree, id) {
    getNodeAge(tree, tree@nodelist[[id]][['prid']])
  }, tree=tree)
  mins <- sapply(ids, getNodeAge, tree=tree)
  data.frame(edge=ids, max=maxs, min=mins, row.names=NULL)
}

# @name getParent
getParent <- function(tree, ids) {
  prids <- getNodesPrid(tree, ids)
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
  pre_1 <- getNodePrid(tree, from)
  pre_2 <- getNodePrid(tree, to)
  parent <- pre_1[which(pre_1 %in% pre_2)[1]]
  path_1 <- c(from ,pre_1[!pre_1 %in% pre_2])
  path_2 <- c(pre_2[!pre_2 %in% pre_1], to)
  c(path_1, parent, path_2)
}

# @name get_Pre
getNodePrid <- function(tree, id) {
  .get <- function(nd, prids) {
    prid <- tree@nodelist[[nd]][['prid']]
    if(!is.null(prid)) {
      prids <- c(prid, .get(prid, prids))
    }
    prids
  }
  .get(id, NULL)
}

getNodesPrid <- function(tree, ids) {
  sapply(ids, getNodePrid, tree=tree, simplify=FALSE)
}

# @name get_Lineage
getNodeLineage <- function(tree, id) {
  prids <- getNodePrid(tree, id)
  lineage <- sapply(prids, function(n) tree@nodelist[[n]][['taxonym']],
                    simplify=FALSE)
  if(length(lineage) > 0) {
    lineage <- c(tree@nodelist[[id]][['taxonym']], lineage)
    lineage <- lineage[length(lineage):1]
    lineage <- unique(lineage)
  } else {
    lineage <- NULL
  }
  lineage
}

getNodesLineage <- function(tree, ids) {
  sapply(ids, getNodeLineage, tree=tree, simplify=FALSE)
}

# @name get_Post
getNodePtid <- function(tree, id) {
  .get <- function(nds, pstids) {
    new_nds <- c()
    for(nd in nds) {
      new_nds <- c(new_nds, tree@nodelist[[nd]][['ptid']])
    }
    pstids <- c(pstids, new_nds)
    if(length(new_nds) > 0) {
      pstids <- .get(nds=new_nds, pstids=pstids)
    }
    pstids
  }
  .get(nds=id, pstids=NULL)
}

getNodesPtid <- function(tree, ids) {
  sapply(ids, getNodePtid, tree=tree)
}

# @name getSubtree
getSubtree <- function(tree, id) {
  pstids <- getNodePtid(tree, id)
  ndlst <- tree@nodelist[c(id, pstids)]
  nd_prdst <- ndlst[[id]][['prdst']]
  ndlst <- lapply(ndlst, function(x) {
    x[['prdst']] <- x[['prdst']] - nd_prdst
    x
  })
  ndlst[[id]][['prid']] <- NULL
  ndlst[[id]][['span']] <- 0
  new_tree <- new('TreeMan', nodelist=ndlst, root=id)
  .update(new_tree)
}
