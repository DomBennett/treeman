# #TODO
setNdTxnym <- function(...) {
  cat('Sorry not yet implemented!\n')
}

setNdsTxnym <- function(...) {
  cat('Sorry not yet implemented!\n')
}

setRoot <- function(...) {
  cat('Sorry not yet implemented!\n')
}

#' @name setPD
#' @title Set the phylogenetic diversity
#' @description Return a tree with the phylogenetic diversity altered.
#' @details Use this function to convert the phylogenetic diversity of a tree. For example,
#' you might want to convert the tree so the sum of all branches is 1. This function will achieve
#' that by modiyfing every branch, while maintaining their relative lengths.
#' @param tree \code{TreeMan} object
#' @param val new phylogenetic diversity
#' @seealso
#' \code{\link{setAge}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' tree <- setPD(tree, val=1)
#' tree <- updateTree(tree)
#' summary(tree)
setPD <- function(tree, val) {
  spns <- getNdsSlt(tree, ids=tree@all, slt_nm="spn")
  spns <- spns/(tree@pd/val)
  tree <- setNdsSpn(tree, ids=tree@all, vals=spns)
  tree@updtd <- FALSE
  tree
}

#' @name setAge
#' @title Set the age of a tree
#' @description Return a tree with the age altered.
#' @details Use this function to change the age of a tree. For example,
#' you might want to convert the tree so that its age equals 1. This function will achieve
#' that by modiyfing every branch, while maintaining their relative lengths.
#' @param tree \code{TreeMan} object
#' @param val new age
#' @seealso
#' \code{\link{setPD}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' tree <- setAge(tree, val=1)
#' tree <- updateTree(tree)
#' summary(tree)
setAge <- function(tree, val) {
  spns <- getNdsSlt(tree, ids=tree@all, slt_nm="spn")
  spns <- spns/(tree@age/val)
  tree <- setNdsSpn(tree, ids=tree@all, vals=spns)
  tree@updtd <- FALSE
  tree
}

#' @name setNdSpn
#' @title Set the branch length of a specific node
#' @description Return a tree with the span of a node altered.
#' @details Takes a tree, a node ID and a new value for the node's preceding branch length (span).
#' Parallelizable.
#' @param tree \code{TreeMan} object
#' @param id id of node whose preceding edge is to be changed
#' @param val new span
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{setNdsSpn}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' tree <- setNdSpn(tree, id='t1', val=100)
#' tree <- updateTree(tree)
setNdSpn <- function(tree, id, val) {
  tree@ndlst[[id]][['spn']] <- val
  tree@updtd <- FALSE
  tree
}

#' @name setNdsSpn
#' @title Set the branch lengths of specific nodes
#' @description Return a tree with the spans of nodes altered.
#' @details Runs \code{setNdSpn} over multiple nodes. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids ids of nodes whose preceding edges are to be changed
#' @param vals new spans
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{setNdSpn}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # make tree taxonomic
#' tree <- setNdsSpn(tree, ids=tree['all'], vals=1)
#' tree <- updateTree(tree)
#' summary(tree)
#' # remove spns by setting all to 0
#' tree <- setNdsSpn(tree, ids=tree['all'], vals=0)
#' tree <- updateTree(tree)
#' summary(tree)
setNdsSpn <- function(tree, ids, vals, ...) {
  .reset <- function(id, spn) {
    ndlst[[id]][['spn']] <- spn
    ndlst[[id]]
  }
  ndlst <- tree@ndlst[ids]
  l_data <- data.frame(id=ids, spn=vals, stringsAsFactors=FALSE)
  ndlst <- plyr::mlply(l_data, .fun=.reset, ...)
  ndlst <- ndlst[1:length(ndlst)]
  tree@ndlst[ids] <- ndlst
  tree@updtd <- FALSE
  tree
}

#' @name setTol
#' @title Set the extinction tolerance
#' @description Return a tree with the tolerance altered.
#' @details Extant tips are determined by how close they are to zero. By default this value
#' is 1e-8. Using this function to change the tolerance will alter the \code{ext} and \code{exc}
#' slots.
#' @param tree \code{TreeMan} object
#' @param tol new tolerance
#' @seealso
#' \code{\link{setNdsSpn}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' tree <- setTol(tree, 10)
#' tree <- updateTree(tree)
#' summary(tree)
setTol <- function(tree, tol) {
  tree@tol <- tol
  tree@updtd <- FALSE
  tree
}

#' @name setNdID
#' @title Set the ID of a node
#' @description Return a tree with the ID of a node altered.
#' @details IDs cannot be changed directly for the \code{TreeMan} class. To change an
#' ID use this function. Warning: all IDs must be unique, avoid spaces in IDs.
#' No need for \code{updateTree()}.
#' @param tree \code{TreeMan} object
#' @param id id to be changed
#' @param val new id
#' @seealso
#' \code{\link{setNdsID}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' tree <- setNdID(tree, 't1', 'heffalump')
setNdID <- function(tree, id, val) {
  setNdsID(tree, id, val)
}

#' @name setNdsID
#' @title Set the IDs of multiple nodes
#' @description Return a tree with the IDs of nodes altered.
#' @details Runs \code{setNdID()} over multiple nodes. Warning: all IDs must be unique,
#' avoid spaces in IDs. Parellizable.
#' No need for \code{updateTree()}.
#' @param tree \code{TreeMan} object
#' @param ids ids to be changed
#' @param vals new ids
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{setNdID}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' new_ids <- paste0('heffalump_', 1:tree['ntips'])
#' tree <- setNdsID(tree, tree['tips'], new_ids)
setNdsID <- function(tree, ids, vals, ...) {
  # internals
  .rplcS4 <- function(slt) {
    if(any(slot(tree, slt) %in% ids)) {
      mtchs <- match(slot(tree, slt), ids)
      return(vals[mtchs])
    } else {
      return(slot(tree, slt))
    }
  }
  .reset <-function(i) {
    .rplc <- function(slt) {
      res <- nd[[slt]]
      mtchs <- match(res, ids)
      res[which(!is.na(mtchs))] <-
        vals[mtchs[!is.na(mtchs)]]
      res
    }
    nd <- tree@ndlst[[i]]
    nd[['id']] <- .rplc("id")
    nd[['ptid']] <- .rplc("ptid")
    nd[['prid']] <- .rplc("prid")
    nd
  }
  l_data <- data.frame(i=1:length(tree@ndlst), stringsAsFactors=FALSE)
  ndlst <- plyr::mlply(l_data, .fun=.reset)
  ndlst <- ndlst[1:length(ndlst)]
  all <- names(tree@ndlst)
  all[match(ids, all)] <- vals
  names(ndlst) <- all
  tree@ndlst <- ndlst
  tree@tips <- .rplcS4('tips')
  tree@nds <- .rplcS4('nds')
  tree@ext <- .rplcS4('ext')
  tree@exc <- .rplcS4('exc')
  tree@root <- .rplcS4('root')
  tree@updtd <- FALSE
  tree
}

#TODO: get these functions working and tested
#' @name setNdOther
#' @title Set a user defined slot
#' @description Return a tree with a user defined slot for node ID.
#' @details A user can specify new slots in a tree. Add a new slot with this function
#' by providing a node ID, a value for the new slot and a unique new slot name. Slot names
#' must not be default \code{TreeMan} names. The new value can be any data type.
#' @param tree \code{TreeMan} object
#' @param id id of the node
#' @param val data for slot
#' @param slt_nm slot name
#' @seealso
#' \code{\link{setNdsOther}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' tree <- setNdOther(tree, 't1', 1, 'binary_val')
#' (getNdSlt(tree, id='t1', slt_nm='binary_val'))
setNdOther <- function(tree, id, val, slt_nm) {
  tree@ndlst[[id]][slt_nm] <- val
  tree@updtd <- FALSE
  tree
}

#' @name setNdsOther
#' @title Set a user defined slot for multiple nodes
#' @description Return a tree with a user defined slot for node IDs.
#' @details Runs \code{setNdOther()} over multiple nodes. Parellizable.
#' @param tree \code{TreeMan} object
#' @param ids id sof the nodes
#' @param vals data for slot
#' @param slt_nm slot name
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{setNdOther}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' vals <- sample(0:1, size=tree['nall'], replace=TRUE)
#' tree <- setNdsOther(tree, tree['all'], vals, 'binary_val')
#' (getNdsSlt(tree, ids=tree['all'], slt_nm='binary_val'))
setNdsOther <- function(tree, ids, vals, slt_nm, ...) {
  .set <- function(id, val) {
    tree@ndlst[[id]][[slt_nm]] <<- val
  }
  l_data <- data.frame(id=ids, val=vals,
                       stringsAsFactors=FALSE)
  plyr::m_ply(.data=l_data, .fun=.set)
  tree@updtd <- FALSE
  tree
}
