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

getNodesSlot <- function(tree, name, ids, ...) {
  .get <- function(i) {
    getNodeSlot(tree, name, ids[i])
  }
  l_data <- data.frame(i=1:length(ids), stringsAsFactors=FALSE)
  res <- mdply(.data=l_data, .fun=.get, ...)
  colnames(res) <- c('id', name)
  res
}

# @name get_Children

getNodeChildren <- function(tree, id) {
  node <- tree@nodelist[[id]]
  node[['children']]
}

getNodesChildren <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- mlply(.data=l_data, .fun=getNodeChildren, tree=tree, ...)
  names(res) <- ids
  res[1:length(res)]
}

# @name get_Age
#TODO: how to effectively handle unrooted trees, age has no meaning
getNodeAge <- function(tree, id) {
  node <- tree@nodelist[[id]]
  age <- tree@age - node[['prdst']]
  age
}

getNodesAge <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- mdply(.data=l_data, .fun=getNodeAge, tree=tree, ...)
  colnames(res) <- c('id', 'age')
  res
}

getSpanAge <- function(tree, id) {
  start <- getNodeAge(tree, tree@nodelist[[id]][['prid']])
  end <- getNodeAge(tree, id)
  data.frame(span=id, start, end)
}

getSpansAge <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- mdply(.data=l_data, .fun=getSpanAge, tree=tree, ...)
  res <- res[ ,colnames(res) != 'id']
  res
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
# recursive, stops whenever prid is NULL
.getNodePrid <- function(tree, id,
                         stopfnc=function(prid){is.null(prid)}) {
  .get <- function(id, prids) {
    prid <- tree@nodelist[[id]][['prid']]
    if(!stopfnc(prid)) {
      prids <- c(prid, .get(prid, prids))
    }
    prids
  }
  .get(id, NULL)
}

getNodePrid <- function(tree, id) {
  .getNodePrid(tree, id)
}

getNodesPrid <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- mlply(.data=l_data, .fun=.getNodePrid, tree=tree, ...)
  names(res) <- ids
  res[1:length(res)]
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

getNodesLineage <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  mlply(.data=l_data, .fun=getNodeLineage, tree=tree, ...)
}

# @name get_Ptid
# reduce dependence on the recursive, by getting prenodes
getNodePtid <- function(tree, id) {
  .get <- function(id, tree) {
    stopfnc <- function(prid) {
      prid %in% pstids
    }
    pstids <<- c(pstids, .getNodePrid(tree, id, stopfnc))
    NULL
  }
  pstids <- id
  l_data <- data.frame(id=tree@nodelist[[id]][['children']],
                       stringsAsFactors=FALSE)
  m_ply(.data=l_data, .fun=.get, tree=tree)
  pstids
}

getNodesPtid <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- mlply(.data=l_data, .fun=getNodePtid, tree=tree, ...)
  names(res) <- ids
  res[1:length(res)]
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
  .updateSlots(new_tree)
}
