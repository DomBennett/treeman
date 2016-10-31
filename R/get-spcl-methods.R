# SINGLE ND
# TODO: bring outgroup, parent and path into terminological line with getNd(s)

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
  # using ndlst guarrantees order
  prids <- .getNdsPridsFrmLst(tree@ndlst, ids=ids, prinds=tree@prinds,
                              parallel=FALSE, progress="none")
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

# SPECIAL

#' @name getTxnyms
#' @title Get IDs for nodes represented txnyms
#' @description Return a list of IDs for any node that contains the given txnyms.
#' @details Returns a list. Txnyms must be spelt correctly.
#' @param tree \code{TreeMan} object
#' @param txnyms vector of taxonomic group names
#' @seealso
#' \code{\link{taxaResolve}}, \code{\link{setTxnyms}}, \code{\link{searchTxnyms}},
#' \code{\link{getNdsLng}}, \code{\link{getNdLng}}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # what ID represents the apes?
#' getTxnyms(mammals, 'Hominoidea')
getTxnyms <- function(tree, txnyms) {
  # get nd id(s) for taxonyms
  .get <- function(id, txnym, ...) {
    for(t in txnyms) {
      if(t %in% txnym) {
        res[[t]] <<- c(res[[t]], id)
      }
    }
  }
  res <- list()
  plyr::m_ply(tree@ndlst, .fun=.get)
  res
}

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
#' apes <- updateTree(apes)
getSubtree <- function(tree, id) {
  if(!id %in% tree@nds) {
    stop('`id` is not an internal node')
  }
  ids <- c(id, getNdPtids(tree, id))
  ndlst <- tree@ndlst[ids]
  ndlst[[id]][['prid']] <- id
  ndlst[[id]][['spn']] <- 0
  new_tree <- new('TreeMan', ndlst=ndlst, root=id,
                  ndmtrx=NULL)
  new_tree
}

# TREE FUNCTIONS

#' @name getTreeAge
#' @title Get age of tree
#' @description Returns age, numeric, of tree
#' @details This can also be achieved with \code{tree['age']} but will
#' only work if the the tree has been updated with \code{updateTree()}.
#' For faster computation, especially within function that perform multiple
#' tree manipulations where the whole tree doesn't need updating, use this function.
#' Parallelizable.
#' @param tree \code{TreeMan} object
#' @param parallel logical, make parallel?
#' @seealso
#' \code{\link{updateTree}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' (getTreeAge(tree))
getTreeAge <- function(tree, parallel=FALSE) {
  tids <- tree@tips
  if(tree@updtd) {
    all_ids <- tree@all
    spns <- .getSltSpns(tree@ndlst)
    res <- .getTreeAgeFrmMtrx(tree@ndmtrx, all_ids, tids, spns, parallel)
  } else {
    res <- .getTreeAgeFrmLst(tree@ndlst, prinds=tree@prinds,
                             tids=tids, parallel)
  }
  res
}