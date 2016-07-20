#' @name getTreeAge
#' @title Get age of tree
#' @description Returns age, numeric, of tree
#' @details This can also be achieved with \code{tree['age']} but will
#' only work if the the tree has been updated with \code{updateTree()}.
#' For faster computation, especially within function that perform multiple
#' tree manipulations where the whole tree doesn't need updating, use this function.
#' @param tree \code{TreeMan} object
#' @seealso
#' \code{\link{updateTree}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' (getTreeAge(tree))
getTreeAge <- function(tree) {
  .getAge(tree@ndlst)
}

#' @name getNdSstr
#' @title Get sister id
#' @description Returns the id of the sister(s) of node id given.
#' @details An error is raised if there is no sister (e.g. for the root).
#'  There can be more than one sister if tree is polytomous.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsSstr}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdSstr(tree, id='t1')
getNdSstr <- function(tree, id) {
  .getSstr(tree@ndlst, id)
}

#' @name getNdsSstr
#' @title Get sister id
#' @description Returns the ids of the sister(s) of nd ids given.
#' @details An error is raised if there is no sister (e.g. for the root).
#'  There can be more than one sister if tree is polytomous. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids nd ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdSstr}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsSstr(tree, ids=tree['tips'])
getNdsSstr <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- plyr::mdply(.data=l_data, .fun=getNdSstr, tree=tree, ...)
  res[ ,2]
}

# #' @name getTxnyms
# #' @title Get node id(s) for txonyms
# #' @description Returns the node ids of nodes with given taxonyms.
# #' @details Returns a \code{list}, parallelizable.
# #' @param tree \code{TreeMan} object
# #' @param txnyms vector of taxonomic names
# #' @param ... \code{plyr} arguments
# #' @seealso
# #' \code{\link{getNdLng}}, 
# #' \code{\link{getNdsLng}}, 
# #' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
# #' @export
# #' @examples
# #' library(treeman)
# #' data(mammals)
# #' homo_ids <- getTxnyms(mammals, txnyms='Homo')
# getTxnyms <- function(tree, txnyms, ...) {
#   # get node id(s) for taxonyms
#   .get <- function(id, txnym, ...) {
#     for(t in txnyms) {
#       if(t %in% txnym) {
#         res[[t]] <<- c(res[[t]], id)
#       }
#     }
#   }
#   res <- vector("list", lenght=length(tree@ndlst))
#   plyr::m_ply(tree@ndlst, .fun=.get)
#   res
# }

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
    kids <- getNdKids(tree, id)
    sum(ids %in% kids)
  }
  prnt <- getPrnt(tree, ids)
  ptids <- tree@ndlst[[prnt]][['ptid']]
  cnts <- sapply(ptids, .cntr)
  outnd <- names(cnts)[which.min(cnts)]
  kids <- getNdKids(tree, outnd)
  if(length(kids) == 0) {
    return(outnd)
  }
  outgroup <- ids[ids %in% kids]
  if(length(outgroup) > 1) {
    return(NULL)
  }
  outgroup
}

#' @name getNdPD
#' @title Get phylogenetic diversity of node
#' @description Return summed value of all descending spns
#' @details Sums the lengths of all descending branches from a node.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsPD}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdPD(tree, id='n1')  # return PD of n1 which in this case is for the whole tree
getNdPD <- function(tree, id) {
  .getPD(tree@ndlst, id)
}

#' @name getNdsPD
#' @title Get phylogenetic diversities of nodes
#' @description Return summed value of all descending spns
#' @details Sums the lengths of all descending branches from a node.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdPD}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsPD(tree, ids=tree['all'])  # return PD of all ids
getNdsPD <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  out <- plyr::mdply(.data=l_data, .fun=getNdPD, tree=tree, ...)
  res <- out[ ,2]
  names(res) <- out[ ,1]
  res
}

#' @name getNdPrdst
#' @title Get pre-distance
#' @description Return root to tip distance (prdst) for \code{id}
#' @details Sums the lengths of all branches from \code{id} to root.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsPrdst}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdPrdst(tree, id='t1')  # return the distance to root from t1
getNdPrdst <- function(tree, id) {
  .getPrdst(tree@ndlst, id)
}

#' @name getNdsPrdst
#' @title Get pre-distances
#' @description Return root to tip distances (prdst) for \code{ids}
#' @details Sums the lengths of all branches from \code{ids} to root.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdPrdst}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsPrdst(tree, ids=tree['tips'])  # return prdsts for all tips
getNdsPrdst <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  out <- plyr::mdply(.data=l_data, .fun=getNdPrdst, tree=tree, ...)
  res <- out[ ,2]
  names(res) <- out[ ,1]
  res
}

#' @name getNdSlt
#' @title Get a node slot
#' @description Returns the value of named slot.
#' @details Returned object depends on name, either character, vector or numeric.
#' Default node slots are: id, spn, prid, ptid and txnym.
#' @param tree \code{TreeMan} object
#' @param slt_nm slot name
#' @param id node id
#' @seealso
#' \code{\link{getNdsSlt}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdSlt(tree, slt_nm='spn', id='t1')  # return span of t1
getNdSlt <- function(tree, slt_nm, id) {
  tree@ndlst[[id]][[slt_nm]]
}

#' @name getNdsSlt
#' @title Get a node slot for multiple nodes
#' @description Returns the value of named slot.
#' @details Returned object depends on name, either character, vector or numeric. Parallelizable.
#' Default node slots are: id, spn, prid, ptid and txnym.
#' @param tree \code{TreeMan} object
#' @param slt_nm slot name
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdSlt}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsSlt(tree, slt_nm='spn', ids=tree['tips'])  # return spans of all tips
getNdsSlt <- function(tree, slt_nm, ids, ...) {
  .get <- function(i) {
    getNdSlt(tree, slt_nm, ids[i])
  }
  l_data <- data.frame(i=1:length(ids), stringsAsFactors=FALSE)
  res <- plyr::mdply(.data=l_data, .fun=.get, ...)
  res[ ,2]
}

#' @name getNdKids
#' @title Get children IDs
#' @description Return the node ids of all tips that descend from node.
#' @details Returns a vector
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsKids}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # everyone descends from root
#' getNdKids(tree, id=tree['root'])
getNdKids <- function(tree, id) {
  names(tree@ndlst)[.getKids(tree@ndlst, id)]
}

#' @name getNdsKids
#' @title Get children IDs for multiple nodes
#' @description Return the node ids of all tips that descend from each node in \code{ids}.
#' @details Returns a list, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdKids}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsKids(tree, id=tree['nds'])
getNdsKids <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- plyr::mlply(.data=l_data, .fun=getNdKids, tree=tree, ...)
  names(res) <- ids
  res[1:length(res)]
}

#' @name getNdAge
#' @title Get age
#' @description Return the age for \code{id}. Requires the known age of the tree to be provided.
#' @details Returns a numeric.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @param tree_age numeric value of known age of tree, tree['age'] if tree is updated
#' @seealso
#' \code{\link{getNdsAge}}, 
#' \code{\link{getSpnAge}}, 
#' \code{\link{getSpnsAge}}, 
#' \code{\link{getPrnt}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # when did apes emerge?
#' # get parent id for all apes
#' prnt_id <- getPrnt(mammals, ids=c('Homo_sapiens', 'Hylobates_concolor'))
#' getNdAge(mammals, id=prnt_id, tree_age=mammals['age'])
getNdAge <- function(tree, id, tree_age) {
  tree_age - .getPrdst(tree@ndlst, id)
}

#' @name getNdsAge
#' @title Get ages for multiple nodes
#' @description Return the age for \code{ids}.
#' @details Returns a vector, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param tree_age numeric value of known age of tree, tree['age'] if tree is updated
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdAge}}, 
#' \code{\link{getSpnAge}}, 
#' \code{\link{getSpnsAge}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsAge(tree, ids=tree['nds'], tree_age=tree['age'])
getNdsAge <- function(tree, ids, tree_age, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  out <- plyr::mdply(.data=l_data, .fun=getNdAge, tree=tree,
                     tree_age=tree_age, ...)
  res <- out[ ,2]
  names(res) <- out[ ,1]
  res
}

#' @name getSpnAge
#' @title Get age range
#' @description Return start and end ages for \code{id} from when it first appears to when it splits
#' @details Returns a dataframe.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @param tree_age numeric value of known age of tree, tree['age'] if tree is updated
#' @seealso
#' \code{\link{getNdAge}}, 
#' \code{\link{getNdsAge}}, 
#' \code{\link{getSpnsAge}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' getSpnAge(mammals, id='Homo_sapiens', tree_age=mammals['age'])
getSpnAge <- function(tree, id, tree_age) {
  start <- getNdAge(tree, tree@ndlst[[id]][['prid']][1], tree_age)
  end <- getNdAge(tree, id, tree_age)
  data.frame(spn=id, start, end)
}

#' @name getSpnsAge
#' @title Get age ranges for multiple nodes
#' @description Return start and end ages for \code{ids} from when they first appear to when they split
#' @details Returns a dataframe, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param tree_age numeric value of known age of tree, tree['age'] if tree is updated
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdAge}}, 
#' \code{\link{getNdsAge}}, 
#' \code{\link{getSpnAge}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # all nodes but root
#' ids <- tree['nds'][tree['nds'] != tree['root']]
#' getSpnsAge(tree, ids=ids, tree_age=tree['age'])
getSpnsAge <- function(tree, ids, tree_age, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- plyr::mdply(.data=l_data, .fun=getSpnAge, tree=tree,
                     tree_age=tree_age, ...)
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
  prids <- getNdsPrids(tree, ids)
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
#' sum(getNdsSlt(mammals, ids=pth, slt_nm='spn'))
getPath <- function(tree, from, to) {
  pre_1 <- c(from, getNdPrids(tree, from))
  pre_2 <- c(to, getNdPrids(tree, to))
  parent <- pre_1[which(pre_1 %in% pre_2)[1]]
  path_1 <- pre_1[!pre_1 %in% pre_2]
  path_2 <- pre_2[!pre_2 %in% pre_1]
  path_2 <- path_2[length(path_2):1]
  c(path_1, parent, path_2)
}

#' @name getNdPrid
#' @title Get pre-nodes
#' @description Return id for directly preceding node
#' @details Returns character.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdPrids}},
#' \code{\link{getNdsPrids}},
#' \code{\link{getNdsPrid}}, 
#' \code{\link{getNdPtid}}, 
#' \code{\link{getNdsPtid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdPrid(tree, id='t1')
getNdPrid <- function(tree, id) {
  tree@ndlst[[id]][['prid']]
}

#' @name getNdsPrid
#' @title Get pre-nodes for multiple nodes
#' @description Return node ids for preceding nodes.
#' @details Returns a list, parallizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdPrids}},
#' \code{\link{getNdsPrids}},
#' \code{\link{getNdPrid}}, 
#' \code{\link{getNdPtid}}, 
#' \code{\link{getNdsPtid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsPrid(tree, ids=tree['tips'])
getNdsPrid <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  out <- plyr::mdply(.data=l_data, .fun=getNdPrid, tree=tree, ...)
  res <- out[ ,2]
  names(res) <- out[ ,1]
  res
}

#' @name getNdPrids
#' @title Get pre-nodes to root
#' @description Return node ids for connecting \code{id} to root.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdPrid}},
#' \code{\link{getNdsPrid}},
#' \code{\link{getNdsPrids}}, 
#' \code{\link{getNdPtid}}, 
#' \code{\link{getNdsPtid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # get all nodes to root
#' getNdPrids(tree, id='t1')
getNdPrids <- function(tree, id) {
  names(tree@ndlst)[.getPrids(tree@ndlst, id)]
}

#' @name getNdsPrids
#' @title Get pre-nodes for multiple nodes
#' @description Return node ids for connecting \code{id} to root.
#' @details Returns a list, parallizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdPrids}},
#' \code{\link{getNdPrid}},
#' \code{\link{getNdsPrid}}, 
#' \code{\link{getNdPtid}}, 
#' \code{\link{getNdsPtid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsPrids(tree, ids=tree['tips'])
getNdsPrids <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  out <- plyr::mlply(.data=l_data, .fun=getNdPrids, tree=tree, ...)
  out <- plyr::mlply(.data=l_data, .fun=getNdPrids, tree=tree)
  names(out) <- attr(out, 'split_labels')[,1]
  res <- out[1:length(out)]
  res
}

#' @name getNdPtid
#' @title Get post-nodes to tips
#' @description Return node ids for connecting \code{id} to kids.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsPtid}}, 
#' \code{\link{getNdPrid}}, 
#' \code{\link{getNdsPrid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # get all nodes from root to tip
#' getNdPtid(tree, id='n1')
# reduce dependence on the recursive, by getting prenodes
# tip ids to id
# TODO: create a ptid and ptids function
getNdPtid <- function(tree, id) {
  names(tree@ndlst)[.getPtids(tree@ndlst, id)]
}

#' @name getNdsPtid
#' @title Get post-nodes to tips for multiple nodes
#' @description Return node ids for connecting \code{ids} to kids.
#' @details Returns a list, parallizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{getNdPtid}}, 
#' \code{\link{getNdPrid}}, 
#' \code{\link{getNdsPrid}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # get all nodes to tip for all nodes
#' getNdsPtid(tree, ids=tree['nds'])
getNdsPtid <- function(tree, ids, ...) {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  out <- plyr::mlply(.data=l_data, .fun=getNdPtid, tree=tree, ...)
  names(out) <- attr(out, 'split_labels')[,1]
  res <- out[1:length(out)]
  res
}

# #' @name getNdLng
# #' @title Get lineage
# #' @description Return unique taxonyms for connecting \code{id} to root.
# #' @details Returns a vector.
# #' @param tree \code{TreeMan} object
# #' @param id node id
# #' @seealso
# #' \code{\link{getNdsLng}}, 
# #' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
# #' @export
# #' @examples
# #' library(treeman)
# #' data(mammals)
# #' # return human lineage
# #' getNdLng(mammals, id='Homo_sapiens')
# getNdLng <- function(tree, id) {
#   .get <- function(txnym, ...) {
#     lng <<- c(txnym, lng)
#   }
#   prids <- c(id, getNdPrids(tree, id))
#   lng <- NULL
#   plyr::m_ply(tree@ndlst[prids], .fun=.get)
#   unique(lng)
# }

# #' @name getNdsLng
# #' @title Get lineage for multiple nodes
# #' @description Return unique taxonyms for connecting \code{ids} to root.
# #' @details Returns a list, parallelizable.
# #' @param tree \code{TreeMan} object
# #' @param ids vector of node ids
# #' @param ... \code{plyr} arguments
# #' @seealso
# #' \code{\link{getNdLng}}, 
# #' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
# #' @export
# #' @examples
# #' library(treeman)
# #' data(mammals)
# #' # return human and gorilla lineages
# #' getNdLng(mammals, id=c('Homo_sapiens', 'Gorilla_gorilla'))
# getNdsLng <- function(tree, ids, ...) {
#   l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
#   out <- plyr::mlply(.data=l_data, .fun=getNdLng, tree=tree, ...)
#   names(out) <- attr(out, 'split_labels')[,1]
#   res <- out[1:length(out)]
#   res
# }

#' @name getSubtree
#' @title Get subtree
#' @description Return tree descending from \code{id}.
#' @details Returns a \code{TreeMan}, parallelizable. \code{id} must be an internal node.
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
  if(!id %in% tree@nds) {
    stop('`id` is not an internal node')
  }
  ids <- c(id, getNdPtid(tree, id))
  ndlst <- tree@ndlst[ids]
  ndlst[[id]][['prid']] <- id
  ndlst[[id]][['spn']] <- 0
  new_tree <- new('TreeMan', ndlst=ndlst, root=id)
  updateTree(new_tree)
}