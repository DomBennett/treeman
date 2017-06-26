#' @name setTxnyms
#' @title Set the txnym slots in a tree
#' @description Return a tree with txnyms added to specified nodes
#' @details Returns a tree. Specify the taxonomic groups for nodes in a tree
#' by providing a vector or list named by node IDs. Takes output from \code{searchTxnyms}.
#' Only letters, numbers and underscores allowed. To remove special characters use regular
#' expressions, e.g. \code{gsub(['a-zA-Z0-9_'], '', txnym)}
#' @param tree \code{TreeMan} object
#' @param txnyms named vector or list
#' @seealso
#' \code{\link{taxaResolve}}, \code{\link{searchTxnyms}},
#' \code{\link{getNdsLng}}, \code{\link{getNdLng}},
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # let's change the txnym for humans
#' # what's its summary before we change anything?
#' summary(mammals[['Homo_sapiens']])
#' # now let's add Hominini
#' new_txnym <- list('Homo_sapiens'=c('Hominini', 'Homo'))
#' mammals <- setTxnyms(mammals, new_txnym)
#' summary(mammals[['Homo_sapiens']])
setTxnyms <- function(tree, txnyms) {
  .add <- function(nid) {
    for(txnym in txnyms[[nid]]) {
      if(grepl('[^a-zA-Z_0-9]', txnym)) {
        stop(paste0('Unsuitable characters in [',
                    txnym, ']'))
      }
    }
    tree@ndlst[[nid]][['txnym']] <<- txnyms[[nid]]
  }
  pull <- names(txnyms) %in% names(tree@ndlst)
  txnyms <- txnyms[pull]
  plyr::m_ply(names(txnyms), .fun=.add)
  tree@wtxnyms <- TRUE
  tree
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
#' summary(tree)
setPD <- function(tree, val) {
  spns <- getNdsSlt(tree, ids=tree@all, slt_nm="spn")
  spns <- spns/(tree@pd/val)
  tree <- setNdsSpn(tree, ids=tree@all, vals=spns)
  tree@pd <- val
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
#' summary(tree)
setAge <- function(tree, val) {
  tree_age <- getAge(tree)
  spns <- getNdsSlt(tree, ids=tree@all, slt_nm="spn")
  spns <- spns/(tree_age)
  tree <- setNdsSpn(tree, ids=tree@all, vals=spns)
  tree
}

#' @name setNdSpn
#' @title Set the branch length of a specific node
#' @description Return a tree with the span of a node altered.
#' @details Takes a tree, a node ID and a new value for the node's preceding branch length (span).
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
#' tree <- updateSlts(tree)
#' summary(tree)
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
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{setNdSpn}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # make tree taxonomic
#' tree <- setNdsSpn(tree, ids=tree['all'], vals=1)
#' summary(tree)
#' # remove spns by setting all to 0
#' tree <- setNdsSpn(tree, ids=tree['all'], vals=0)
#' summary(tree)
setNdsSpn <- function(tree, ids, vals, parallel=FALSE, progress="none") {
  .reset <- function(id, spn) {
    ndlst[[id]][['spn']] <- spn
    ndlst[[id]]
  }
  ndlst <- tree@ndlst[ids]
  l_data <- data.frame(id=ids, spn=vals, stringsAsFactors=FALSE)
  ndlst <- plyr::mlply(l_data, .fun=.reset, .parallel=parallel,
                       .progress=progress)
  ndlst <- ndlst[1:length(ndlst)]
  tree@ndlst[ids] <- ndlst
  tree <- updateSlts(tree)
  tree
}

#' @name setNdID
#' @title Set the ID of a node
#' @description Return a tree with the ID of a node altered.
#' @details IDs cannot be changed directly for the \code{TreeMan} class. To change an
#' ID use this function. Warning: all IDs must be unique, avoid spaces in IDs and only
#' use letters, numbers and underscores.
#' Use \code{\link{updateSlts}} after running.
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
#' tree <- updateSlts(tree)
setNdID <- function(tree, id, val) {
  tree@updtd <- FALSE
  setNdsID(tree, id, val)
}

#' @name setNdsID
#' @title Set the IDs of multiple nodes
#' @description Return a tree with the IDs of nodes altered.
#' @details Runs \code{setNdID()} over multiple nodes. Warning: all IDs must be unique,
#' avoid spaces in IDs, only use numbers, letters and underscores. Parellizable.
#' @param tree \code{TreeMan} object
#' @param ids ids to be changed
#' @param vals new ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{setNdID}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' new_ids <- paste0('heffalump_', 1:tree['ntips'])
#' tree <- setNdsID(tree, tree['tips'], new_ids)
#' summary(tree)
setNdsID <- function(tree, ids, vals, parallel=FALSE, progress="none") {
  # internals
  .testSpcls <- function(id) {
    if(grepl('[^a-zA-Z_0-9]', id)) {
      stop(paste0('Unsuitable characters in [', id, ']'))
    }
    NULL
  }
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
  sapply(vals, .testSpcls)
  l_data <- data.frame(i=1:length(tree@ndlst), stringsAsFactors=FALSE)
  ndlst <- plyr::mlply(l_data, .fun=.reset, .parallel=parallel, .progress=progress)
  ndlst <- ndlst[1:length(ndlst)]
  all <- names(tree@ndlst)
  all[match(ids, all)] <- vals
  names(ndlst) <- all
  tree@ndlst <- ndlst
  tree@tips <- .rplcS4('tips')
  tree@nds <- .rplcS4('nds')
  tree@root <- .rplcS4('root')
  tree <- updateSlts(tree)
  tree
}

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
#' tree <- updateSlts(tree)
#' (getNdSlt(tree, id='t1', slt_nm='binary_val'))
setNdOther <- function(tree, id, val, slt_nm) {
  tree@ndlst[[id]][slt_nm] <- val
  tree@updtd <- FALSE
  if(!slt_nm %in% tree@othr_slt_nms) {
    tree@othr_slt_nms <- c(tree@othr_slt_nms, slt_nm)
  }
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
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{setNdOther}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # e.g. confidences for nodes
#' vals <- runif(min=0, max=1, n=tree['nall'])
#' tree <- setNdsOther(tree, tree['all'], vals, 'confidence')
#' tree <- updateSlts(tree)
#' summary(tree)
#' (getNdsSlt(tree, ids=tree['all'], slt_nm='confidence'))
setNdsOther <- function(tree, ids, vals, slt_nm, parallel=FALSE, progress="none") {
  .set <- function(id, val) {
    tree@ndlst[[id]][[slt_nm]] <<- val
  }
  l_data <- data.frame(id=ids, val=vals,
                       stringsAsFactors=FALSE)
  plyr::m_ply(.data=l_data, .fun=.set, .parallel=parallel, .progress=progress)
  tree@updtd <- FALSE
  if(!slt_nm %in% tree@othr_slt_nms) {
    tree@othr_slt_nms <- c(tree@othr_slt_nms, slt_nm)
  }
  tree
}

#' @name rmOtherSlt
#' @title Remove a user-defined slot
#' @description Returns a tree with a user-defined tree slot removed.
#' @details A user can specify a new slot using the \code{setNdSlt()} function
#' or upon reading a tree. This can be removed using this function by specifying
#' the name of the slot to be removed.
#' @param tree \code{TreeMan} object
#' @param slt_nm name of slot to be removed
#' @seealso
#' \code{\link{setNdOther}}, \code{\link{setNdsOther}},
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' vals <- runif(min=0, max=1, n=tree['nall'])
#' tree <- setNdsOther(tree, tree['all'], vals, 'confidence')
#' tree <- updateSlts(tree)
#' summary(tree)
#' tree <- rmOtherSlt(tree, 'confidence')
#' tree <- updateSlts(tree)
#' summary(tree)
rmOtherSlt <- function(tree, slt_nm) {
  .set <- function(id) {
    tree@ndlst[[id]][[slt_nm]] <<- NULL
  }
  l_data <- data.frame(id=tree@all, stringsAsFactors=FALSE)
  plyr::m_ply(.data=l_data, .fun=.set)
  tree@updtd <- FALSE
  tree@othr_slt_nms <- tree@othr_slt_nms[tree@othr_slt_nms != slt_nm]
  tree
}
