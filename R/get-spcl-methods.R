# bipartitions
#' @name getBiprts
#' @title Get the sets of labels for each bipartition in tree
#' @description Returns a list of tip IDs for each branch in the tree. Options
#' allow the user to act as if the root is not present and to use a universal
#' code for comparing between trees.
#' @details Setting \code{root} to FALSE will ignore the bipartitions created by
#' the root. Setting \code{universal} to TRUE will return a vector of 0s and 1s,
#' not a list of tips. These codes will always begin with 1, and will allow for
#' the comparison of splits between trees as they do not have "chiralty", so to
#' speak.
#' @param tree \code{TreeMan} object
#' @param tips vector of tips IDs to use for bipartitions
#' @param root Include the root for the bipartitions? Default TRUE.
#' @param universal Create a code for comparing between trees
#' @seealso
#' \code{\link{calcDstRF}}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # get all of the tip IDs for each branch in the rooted tree
#' (getBiprts(tree))
#' # ignore the root and get bipartitions for unrooted tree
#' (getBiprts(tree, root = FALSE))
#' # use the universal code for comparing splits between trees
#' (getBiprts(tree, root = FALSE, universal = TRUE))
getBiprts <- function(tree, tips = tree@tips, root = TRUE, universal = FALSE) {
  kids <- getNdsKids(tree = tree, ids = tree@nds)
  res <- lapply(X = kids, FUN = function(x) tips %in% x)
  if (!root) {
    n <- vapply(X = res, FUN = sum, FUN.VALUE = integer(1))
    # drop splits consisting of all tips or just 1
    res <- res[n < (length(tips) - 1) & n > 1]
  }
  if (universal) {
    res <- unname(vapply(X = res, FUN = function(x) {
      if (!x[[1]]) {
        x <- !x
      }
      paste0(as.integer(x), collapse = '')
    }, FUN.VALUE = character(1)))
    res <- unique(res)
  } else {
    res <- lapply(X = res, FUN = function(x) {
      list(tree@tips[x], tree@tips[!x])
    })
  }
  res
}

# ULTRAMETRIC
#' @name isUltrmtrc
#' @title Is tree ultrametric?
#' @description Return TRUE if all tips end at 0, else FALSE.
#' @details Returns a boolean. This function works in the background
#' for the \code{['ultr']} slot in a \code{TreeMan} object.
#' @param tree \code{TreeMan} object
#' @param tol zero tolerance
#' @seealso
#' \code{\link{getLvng}}, \code{\link{getDcsd}}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' (isUltrmtrc(tree))
isUltrmtrc <- function(tree, tol=1e-8) {
  dead <- .livingOrDesceased(tree, tol=tol, bool=FALSE)
  if(length(dead) > 0) {
    return(FALSE)
  }
  TRUE
}

# EXTINCT/EXTANT
.livingOrDesceased <- function(tree, tol=1e-8, bool) {
  if(!is.null(tree@ndmtrx)) {
    spns <- getNdsSlt(tree, 'spn', names(tree@ndlst))
    tip_prdsts <- .getNdsPrdstsFrmMtrx(tree@ndmtrx, tree@all,
                                       tree@tips, spns,
                                       parallel=FALSE,
                                       progress="none")
  } else {
    tip_prdsts <- .getNdsPrdstsFrmLst(tree@ndlst, tree@tips,
                                      tree@prinds, parallel=FALSE,
                                      progress='none')
  }
  age <- max(tip_prdsts)
  extant_is <- (age - tip_prdsts) <= tol
  living <- names(extant_is)[extant_is]
  deceased <- tree@tips[!tree@tips %in% living]
  if(bool) {
    return(living)
  }
  deceased
}

#' @name getDcsd
#' @title Get extinct tips from a tree
#' @description Return all extinct tip \code{ID}s.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param tol zero tolerance
#' @seealso
#' \code{\link{getLvng}}, \code{\link{isUltrmtrc}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' (getDcsd(tree))
getDcsd <- function(tree, tol=1e-8) {
  .livingOrDesceased(tree=tree, tol=tol, bool=FALSE)
}

#' @name getLvng
#' @title Get extant tips from a tree
#' @description Return all extant tip \code{ID}s.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param tol zero tolerance
#' @seealso
#' \code{\link{getDcsd}}, \code{\link{isUltrmtrc}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' (getLvng(tree))
getLvng <- function(tree, tol=1e-8) {
  .livingOrDesceased(tree=tree, tol=tol, bool=TRUE)
}

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
  cnts <- vapply(ptids, .cntr, integer(1))
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

#' @name getUnqNds
#' @title Get unique nodes represented by tips
#' @description Return a list of IDs for any node that are represented by tip IDs given.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param tids vector of tip IDs
#' @seealso
#' \code{\link{getCnnctdNds}}, \code{\link{calcFrPrp}},
#' \code{\link{calcPhyDv}}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' unqnds <- getUnqNds(tree, c('t1', 't2'))
getUnqNds <- function(tree, tids) {
  rmng <- tree@tips[!tree@tips %in% tids]
  ignr <- c(unique(unlist(getNdsPrids(tree, tids))), tids)
  rmng <- c(unique(unlist(getNdsPrids(tree, rmng))), rmng)
  ignr[!ignr %in% rmng]
}

#' @name getCnnctdNds
#' @title Get all nodes connected by given tips
#' @description Return a vector of IDs of all nodes that are connected to tip IDs given.
#' @details Returns a vector. This function is the basis for \code{calcPhyDv()}, it determines
#' the unique set of nodes connected for a set of tips.
#' @param tree \code{TreeMan} object
#' @param tids vector of tip IDs
#' @seealso
#' \code{\link{getUnqNds}}, \code{\link{calcFrPrp}},
#' \code{\link{calcPhyDv}}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' cnntdnds <- getCnnctdNds(tree, c('t1', 't2'))
getCnnctdNds <- function(tree, tids) {
  prids <- c(unlist(getNdsPrids(tree, tids)),
             tids)
  counts <- table(prids)
  names(counts)[counts < length(tids)]
}

#' @name getNdsFrmTxnyms
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
#' getNdsFrmTxnyms(mammals, 'Hominoidea')
getNdsFrmTxnyms <- function(tree, txnyms) {
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
#' \code{\link{getPrnt}}, \code{\link{addClade}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # get tree of apes
#' ape_id <- getPrnt(mammals, ids=c('Homo_sapiens', 'Hylobates_concolor'))
#' apes <- getSubtree(mammals, id=ape_id)
#' summary(apes)
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
  new_tree <- pstMnp(new_tree)
  new_tree <- updateSlts(new_tree)
  new_tree
}

# TREE FUNCTIONS

#' @name getAge
#' @title Get age of tree
#' @description Returns age, numeric, of tree
#' @details Calculates the age of a tree, determined as the maximum tip to root
#' distance.
#' @param tree \code{TreeMan} object
#' @param parallel logical, make parallel?
#' @seealso
#' \code{\link{updateSlts}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' (getAge(tree))
getAge <- function(tree, parallel=FALSE) {
  tids <- tree@tips
  if(!is.null(tree@ndmtrx)) {
    all_ids <- tree@all
    spns <- .getSltSpns(tree@ndlst)
    res <- .getTreeAgeFrmMtrx(tree@ndmtrx, all_ids, tids, spns, parallel)
  } else {
    res <- .getTreeAgeFrmLst(tree@ndlst, prinds=tree@prinds,
                             tids=tids, parallel)
  }
  res
}