# TODO: sort documentation

getOutgroup <- function(tree, ids) {
  .cntr <- function(id) {
    kids <- tree@nodelist[[id]][['kids']]
    sum(ids %in% kids)
  }
  prnt <- getParent(tree, ids)
  ptids <- tree@nodelist[[prnt]][['ptid']]
  cnts <- sapply(ptids, .cntr)
  outnd <- names(cnts)[which.min(cnts)]
  kids <- tree@nodelist[[outnd]][['kids']]
  if(length(kids) == 0) {
    return(outnd)
  }
  outgroup <- ids[ids %in% kids]
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
  res[ ,2]
}

# @name get_kids

getNodeKids <- function(tree, id) {
  node <- tree@nodelist[[id]]
  node[['kids']]
}

getNodesKids <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- mlply(.data=l_data, .fun=getNodeKids, tree=tree, ...)
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
  res[ ,2]
}

getSpanAge <- function(tree, id) {
  start <- getNodeAge(tree, tree@nodelist[[id]][['prid']][1])
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
  pre_1 <- c(from, getNodePrid(tree, from))
  pre_2 <- c(to, getNodePrid(tree, to))
  parent <- pre_1[which(pre_1 %in% pre_2)[1]]
  path_1 <- pre_1[!pre_1 %in% pre_2]
  path_2 <- pre_2[!pre_2 %in% pre_1]
  path_2 <- path_2[length(path_2):1]
  c(path_1, parent, path_2)
}

# @name get_Prid
# return from id to stop_id(s) (usually root)
getNodePrid <- function(tree, id, stop_id=tree@root) {
  tree@nodelist[[id]][['prid']]
}

getNodesPrid <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- mlply(.data=l_data, .fun=getNodePrid, tree=tree, ...)
  names(res) <- ids
  res[1:length(res)]
}

# @name get_Ptid
# reduce dependence on the recursive, by getting prenodes
# tip ids to id
getNodePtid <- function(tree, id) {
  .get <- function(id) {
    tmp <- getNodePrid(tree, id, stop_id=pstids)
    tmp <- tmp[-length(tmp)]
    pstids <<- c(tmp, pstids)
    NULL
  }
  pstids <- id
  l_data <- data.frame(id=tree@nodelist[[id]][['children']],
                       stringsAsFactors=FALSE)
  m_ply(.data=l_data, .fun=.get)
  pstids
}

getNodesPtid <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- mlply(.data=l_data, .fun=getNodePtid, tree=tree, ...)
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

# @name getSubtree
getSubtree <- function(tree, id) {
  .prdst <- function(nd) {
    nd[['prdst']] <- nd[['prdst']] - nd_prdst
    nd
  }
  pstids <- getNodePtid(tree, id)
  ndlst <- tree@nodelist[pstids]
  nd_prdst <- ndlst[[id]][['prdst']]
  ndlst <- llply(.data=ndlst, .fun=.prdst)
  ndlst[[id]][['prid']] <- NULL
  ndlst[[id]][['span']] <- 0
  new_tree <- new('TreeMan', nodelist=ndlst, root=id)
  .updateSlots(new_tree)
}
