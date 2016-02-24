#' @name getTxnyms
#' @title Get node id(s) for txonyms
#' @description Returns the node ids of nodes with given taxonyms.
#' @details Returns a \code{list}, parallelizable.
#' @param tree \code{TreeMan} object
#' @param txnyms vector of taxonomic names
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNodeLineage}}, 
#' \code{\link{getNodesLineage}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' homo_ids <- getTxnyms(mammals, txnyms='Homo')
getTxnyms <- function(tree, txnyms, ...) {
  # get node id(s) for taxonyms
  .get <- function(id, txnym, ...) {
    for(t in txnyms) {
      if(t %in% txnym) {
        res[[t]] <<- c(res[[t]], id)
      }
    }
  }
  res <- list()
  m_ply(tree@nodelist, .fun=.get)
  res
}

#' @name getOtgrp
#' @title Get outgroup
#' @description Return the outgroup based on a tree and a vector of IDs.
#' @details Returns a id, character. If there are multiple possible outgroups, returns NULL.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @seealso
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # orangutan is an outgroup wrt humans and chimps
#' getOtgrp(mammals, ids=c('Homo_sapiens', 'Pan_troglodytes', 'Pongo_pygmaeus'))
getOtgrp <- function(tree, ids) {
  .cntr <- function(id) {
    kids <- tree@nodelist[[id]][['kids']]
    sum(ids %in% kids)
  }
  prnt <- getPrnt(tree, ids)
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

#' @name getNodeSlot
#' @title Get a node slot
#' @description Returns the value of named slot.
#' @details Returned object depends on name, either character, vector or numeric.
#' @param tree \code{TreeMan} object
#' @param name slot name
#' @param id node id
#' @seealso
#' \code{\link{getNodesSlot}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNodeSlot(tree, name='span', id='t1')  # return span of t1
getNodeSlot <- function(tree, name, id) {
  tree@nodelist[[id]][[name]]
}

#' @name getNodesSlot
#' @title Get a node slot for multiple nodes
#' @description Returns the value of named slot.
#' @details Returned object depends on name, either character, vector or numeric. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param name slot name
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNodeSlot}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNodesSlot(tree, name='span', ids=tree['tips'])  # return spans of all tips
getNodesSlot <- function(tree, name, ids, ...) {
  .get <- function(i) {
    getNodeSlot(tree, name, ids[i])
  }
  l_data <- data.frame(i=1:length(ids), stringsAsFactors=FALSE)
  res <- mdply(.data=l_data, .fun=.get, ...)
  res[ ,2]
}

#' @name getNodeKids
#' @title Get children IDs
#' @description Return the node ids of all tips that descend from node.
#' @details Returns a vector
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNodesKid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # everyone descends from root
#' getNodeKids(tree, id=tree['root'])

getNodeKids <- function(tree, id) {
  node <- tree@nodelist[[id]]
  node[['kids']]
}

#' @name getNodesKids
#' @title Get children IDs for multiple nodes
#' @description Return the node ids of all tips that descend from each node in \code{ids}.
#' @details Returns a list, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNodeKid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNodesKids(tree, id=tree['nodes'])

getNodesKids <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- mlply(.data=l_data, .fun=getNodeKids, tree=tree, ...)
  names(res) <- ids
  res[1:length(res)]
}

#' @name getNodeAge
#' @title Get age
#' @description Return the root to tip distance for \code{id}.
#' @details Returns a numeric.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNodesAge}}, 
#' \code{\link{getSpanAge}}, 
#' \code{\link{getSpansAge}}, 
#' \code{\link{getPrnt}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # when did apes emerge?
#' prnt_id <- getPrnt(mammals, ids=c('Homo_sapiens', 'Hylobates_concolor'))  # get parent id for all apes
#' getNodeAge(mammals, id=prnt_id)
#TODO: how to effectively handle unrooted trees, age has no meaning
getNodeAge <- function(tree, id) {
  node <- tree@nodelist[[id]]
  age <- tree@age - node[['prdst']]
  age
}

#' @name getNodesAge
#' @title Get ages for multiple nodes
#' @description Return the root to tip distances for \code{ids}.
#' @details Returns a vector, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNodeAge}}, 
#' \code{\link{getSpanAge}}, 
#' \code{\link{getSpansAge}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNodesAge(tree, ids=tree['nodes'])
getNodesAge <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- mdply(.data=l_data, .fun=getNodeAge, tree=tree, ...)
  res[ ,2]
}

#' @name getSpanAge
#' @title Get age range
#' @description Return start and end root to tip distances for \code{id}.
#' @details Returns a dataframe.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNodeAge}}, 
#' \code{\link{getNodesAge}}, 
#' \code{\link{getSpansAge}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' getSpanAge(mammals, id='Homo_sapiens')
getSpanAge <- function(tree, id) {
  start <- getNodeAge(tree, tree@nodelist[[id]][['prid']][1])
  end <- getNodeAge(tree, id)
  data.frame(span=id, start, end)
}

#' @name getSpansAge
#' @title Get age ranges for multiple nodes
#' @description Return start and end root to tip distances for \code{ids}.
#' @details Returns a dataframe, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNodeAge}}, 
#' \code{\link{getNodesAge}}, 
#' \code{\link{getSpanAge}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # all nodes but root
#' ids <- tree['nodes'][tree['nodes'] != tree['root']]
#' getSpansAge(tree, ids=ids)
getSpansAge <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- mdply(.data=l_data, .fun=getSpanAge, tree=tree, ...)
  res <- res[ ,colnames(res) != 'id']
  res
}

#' @name getPrnt
#' @title Get parent
#' @description Return parental (most recent common ancestor) node id for \code{ids}.
#' @details Returns a character.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @seealso
#' \code{\link{getSubtree}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # choosing ids from the two main branches of apes allows to find the parent for all apes
#' ape_id <- getPrnt(mammals, ids=c('Homo_sapiens', 'Hylobates_concolor'))
getPrnt <- function(tree, ids) {
  prids <- getNodesPrid(tree, ids)
  rf <- prids[[1]]
  mn_rnk <- 0
  for(n in prids[-1]) {
    rnk <- min(match(n, rf), na.rm=TRUE)
    if(rnk > mn_rnk) mn_rnk <- rnk
  }
  rf[mn_rnk]
}

#' @name getPath
#' @title Get path between nodes
#' @description Return node ids for connecting \code{from} to \code{to}.
#' @details Returns a vector, first id is \code{from} to \code{to}.
#' @param tree \code{TreeMan} object
#' @param from starting node id
#' @param to ending node id
#' @seealso
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # what's the phylogenetic distance from humans to gorillas?
#' ape_id <- getPrnt(mammals, ids=c('Homo_sapiens', 'Hylobates_concolor'))
#' pth <- getPath(mammals, from='Homo_sapiens', to='Gorilla_gorilla')
#' sum(getNodesSlot(mammals, ids=pth, name='span'))
getPath <- function(tree, from, to) {
  pre_1 <- c(from, getNodePrid(tree, from))
  pre_2 <- c(to, getNodePrid(tree, to))
  parent <- pre_1[which(pre_1 %in% pre_2)[1]]
  path_1 <- pre_1[!pre_1 %in% pre_2]
  path_2 <- pre_2[!pre_2 %in% pre_1]
  path_2 <- path_2[length(path_2):1]
  c(path_1, parent, path_2)
}

#' @name getNodePrid
#' @title Get pre-nodes to root
#' @description Return node ids for connecting \code{id} to root.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNodesPrid}}, 
#' \code{\link{getNodePtid}}, 
#' \code{\link{getNodesPtid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # get all nodes to root
#' getNodePrid(tree, id='t1')
getNodePrid <- function(tree, id) {
  tree@nodelist[[id]][['prid']]
}

#' @name getNodesPrid
#' @title Get pre-nodes to root for multiple nodes
#' @description Return node ids for connecting \code{ids} to root.
#' @details Returns a list, parallizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNodePrid}}, 
#' \code{\link{getNodePtid}}, 
#' \code{\link{getNodesPtid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # get all nodes to root
#' getNodesPrid(tree, ids=tree['tips'])
getNodesPrid <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- mlply(.data=l_data, .fun=getNodePrid, tree=tree, ...)
  names(res) <- ids
  res[1:length(res)]
}

#' @name getNodePtid
#' @title Get post-nodes to tips
#' @description Return node ids for connecting \code{id} to kids.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNodesPtid}}, 
#' \code{\link{getNodePrid}}, 
#' \code{\link{getNodesPrid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # get all nodes from root to tip
#' getNodePtid(tree, id='n1')
# reduce dependence on the recursive, by getting prenodes
# tip ids to id
getNodePtid <- function(tree, id) {
  .get <- function(id) {
    tmp <- c(id, getNodePrid(tree, id))
    index <- seq(from=1, to=(which(tmp %in% pstids)[1]-1), by=1)
    pstids <<- c(tmp[index], pstids)
    NULL
  }
  pstids <- id
  l_data <- data.frame(id=tree@nodelist[[id]][['kids']],
                       stringsAsFactors=FALSE)
  m_ply(.data=l_data, .fun=.get)
  pstids
}

#' @name getNodesPtid
#' @title Get post-nodes to tips for multiple nodes
#' @description Return node ids for connecting \code{ids} to kids.
#' @details Returns a list, parallizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNodePtid}}, 
#' \code{\link{getNodePrid}}, 
#' \code{\link{getNodesPrid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # get all nodes to tip for all nodes
#' getNodesPtid(tree, ids=tree['nodes'])
getNodesPtid <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- mlply(.data=l_data, .fun=getNodePtid, tree=tree, ...)
  names(res) <- ids
  res[1:length(res)]
}

#' @name getNodeLineage
#' @title Get lineage
#' @description Return unique taxonyms for connecting \code{id} to root.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNodesLineage}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # return human lineage
#' getNodeLineage(mammals, id='Homo_sapiens')
getNodeLineage <- function(tree, id) {
  .get <- function(txnym, ...) {
    lng <<- c(txnym, lng)
  }
  prids <- c(id, getNodePrid(tree, id))
  lng <- NULL
  m_ply(tree@nodelist[prids], .fun=.get)
  unique(lng)
}

#' @name getNodesLineage
#' @title Get lineage for multiple nodes
#' @description Return unique taxonyms for connecting \code{ids} to root.
#' @details Returns a list, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNodeLineage}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # return human and gorilla lineages
#' getNodeLineage(mammals, id=c('Homo_sapiens', 'Gorilla_gorilla'))
getNodesLineage <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  mlply(.data=l_data, .fun=getNodeLineage, tree=tree, ...)
}

#' @name getSubtree
#' @title Get subtree
#' @description Return tree descending from \code{id}.
#' @details Returns a \code{TreeMan}, parallelizable.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getPrnt}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # get tree of apes
#' ape_id <- getPrnt(mammals, ids=c('Homo_sapiens', 'Hylobates_concolor'))
#' apes <- getSubtree(mammals, id=ape_id)
getSubtree <- function(tree, id) {
  .prdst <- function(nd) {
    nd[['prdst']] <- nd[['prdst']] - nd_prdst
    nd[['prid']] <- nd[['prid']][nd[['prid']] %in% nids]
    nd
  }
  pstids <- getNodePtid(tree, id)
  ndlst <- tree@nodelist[pstids]
  nids <- names(ndlst)
  nd_prdst <- ndlst[[id]][['prdst']]
  ndlst <- llply(.data=ndlst, .fun=.prdst)
  ndlst[[id]][['prid']] <- NULL
  ndlst[[id]][['span']] <- 0
  new_tree <- new('TreeMan', nodelist=ndlst, root=id)
  .updateTreeSlots(new_tree)
}