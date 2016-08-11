#' @name getNdsSstr
#' @title Get sister id
#' @description Returns the ids of the sister(s) of nd ids given.
#' @details An error is raised if there is no sister (e.g. for the root).
#'  There can be more than one sister if tree is polytomous. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids nd ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdSstr}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsSstr(tree, ids=tree['tips'])
getNdsSstr <- function(tree, ids, parallel=FALSE, progress="none") {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  res <- plyr::mdply(.data=l_data, .fun=.getNdSstrFrmLst, ndlst=tree@ndlst,
                     .parallel=parallel, .progress=progress)
  res[ ,2]
}

#' @name getNdsPD
#' @title Get phylogenetic diversities of nodes
#' @description Return summed value of all descending spns
#' @details Sums the lengths of all descending branches from a node.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param parallel logical, make parallel? #' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdPD}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsPD(tree, ids=tree['all'])  # return PD of all ids
getNdsPD <- function(tree, ids, parallel=FALSE, progress="none") {
  if(tree@updtd & length(ids) > 1) {
    res <- .getNdsPDFrmMtrx(tree@ndsmtrx, ids, parallel, progress)
  } else {
    res <- .getNdsPDFrmLst(tree@ndslst, ids, parallel, progress)
  }
  res
}

#' @name getNdsPrdst
#' @title Get pre-distances
#' @description Return root to tip distances (prdst) for \code{ids}
#' @details Sums the lengths of all branches from \code{ids} to root.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdPrdst}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsPrdst(tree, ids=tree['tips'])  # return prdsts for all tips
getNdsPrdst <- function(tree, ids, parallel=FALSE, progress="none") {
  if(tree@updtd & length(ids) > 1) {
    res <- .getNdsPrdstFrmMtrx(tree@ndsmtrx, ids, parallel, progress)
  } else {
    res <- .getNdsPrdstFrmLst(tree@ndslst, ids, parallel, progress)
  }
  res
}

#' @name getNdsSlt
#' @title Get a node slot for multiple nodes
#' @description Returns the value of named slot.
#' @details Returned object depends on name, either character, vector or numeric. Parallelizable.
#' Default node slots are: id, spn, prid, ptid and txnym.
#' @param tree \code{TreeMan} object
#' @param slt_nm slot name
#' @param ids vector of node ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdSlt}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsSlt(tree, slt_nm='spn', ids=tree['tips'])  # return spans of all tips
getNdsSlt <- function(tree, slt_nm, ids, parallel=FALSE, progress="none") {
  .get <- function(i) {
    getNdSlt(tree, slt_nm, ids[i])
  }
  l_data <- data.frame(i=1:length(ids), stringsAsFactors=FALSE)
  res <- plyr::mdply(.data=l_data, .fun=.get, .parallel=parallel,
                     .progress=progress)
  res[ ,2]
}

#' @name getNdsKids
#' @title Get children IDs for multiple nodes
#' @description Return the node ids of all tips that descend from each node in \code{ids}.
#' @details Returns a list, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdKids}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsKids(tree, id=tree['nds'])
getNdsKids <- function(tree, ids, parallel=FALSE,
                       progress="none") {
  if(tree@updtd & length(ids) > 1) {
    res <- .getNdsKidsFrmMtrx(tree@ndmtrx, tree@all,
                              ids, tree@tips, parallel, progress)
  } else {
    res <- .getNdsKidsFrmLst(tree@ndlst, ids,
                             parallel, progress)
  }
  res
}

#' @name getNdsAge
#' @title Get ages for multiple nodes
#' @description Return the age for \code{ids}.
#' @details Returns a vector, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param tree_age numeric value of known age of tree, tree['age'] if tree is updated
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
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
getNdsAge <- function(tree, ids, tree_age,
                      parallel=FALSE,
                      progress="none") {
  if(tree@updtd & length(ids) > 1) {
    spns <- .getSltSpns(tree@ndlst, parallel)
    res <- .getNdsPrdstsFrmMtrx(tree@ndmtrx, tree@all,
                               ids, spns, parallel, progress)
    res <- tree_age - res
  } else {
    res <- .getNdsPrdstsFrmLst(tree@ndlst, ids,
                              parallel, progress)
    res <- tree_age - res
  }
  res
}

#' @name getSpnsAge
#' @title Get age ranges for multiple nodes
#' @description Return start and end ages for \code{ids} from when they first appear to when they split
#' @details Returns a dataframe, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param tree_age numeric value of known age of tree, tree['age'] if tree is updated
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
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
getSpnsAge <- function(tree, ids, tree_age,
                       parallel=FALSE, progress="none") {
  spns <- .getSltSpns(tree@ndlst, parallel)
  if(tree@updtd & length(ids) > 1) {
    end <- .getNdsPrdstsFrmMtrx(tree@ndmtrx, tree@all,
                                ids, spns, parallel, progress)
    
  } else {
    end <- .getNdsPrdstsFrmLst(tree@ndlst, ids,
                               parallel, progress)
  }
  start <- end - spns
  end <- tree_age - end
  start <- tree_age - start
  data.frame(spn=ids, start, end, row.names=NULL)
}

#' @name getNdsPrids
#' @title Get pre-nodes for multiple nodes
#' @description Return node ids for connecting \code{id} to root.
#' @details Returns a list, parallizable. The function will work faster
#' if \code{ordrd} is FALSE.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ordrd logical, ensure returned prids are ordered ID to root
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdPrids}},
#' \code{\link{getNdPtids}}, 
#' \code{\link{getNdsPtids}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdsPrids(tree, ids=tree['tips'])
getNdsPrids <- function(tree, ids, ordrd=FALSE,
                        parallel=FALSE, progress="none") {
  if(tree@updtd & length(ids) > 1 & !ordrd) {
    res <- .getNdsPridsFrmMtrx(tree@ndmtrx, tree@all,
                               ids, parallel, progress)
  } else {
    res <- .getNdsPridsFrmLst(tree@ndlst, ids,
                              parallel, progress)
  }
  res
}

#' @name getNdsPtids
#' @title Get post-nodes to tips for multiple nodes
#' @description Return node ids for connecting \code{ids} to kids.
#' @details Returns a list, parallizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdPtids}}, 
#' \code{\link{getNdPrids}}, 
#' \code{\link{getNdsPrids}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # get all nodes to tip for all nodes
#' getNdsPtids(tree, ids=tree['nds'])
getNdsPtids <- function(tree, ids, parallel=FALSE, progress="none") {
  if(tree@updtd & length(ids) > 1) {
    res <- .getNdsPtidsFrmMtrx(tree@ndmtrx, tree@all,
                               ids, parallel, progress)
  } else {
    res <- .getNdsPtidsFrmLst(tree@ndlst, ids,
                              parallel, progress)
  }
  res
}

# FOR RESURRECTION
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
# #' @param parallel logical, make parallel? #' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
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

# #' @name getTxnyms
# #' @title Get node id(s) for txonyms
# #' @description Returns the node ids of nodes with given taxonyms.
# #' @details Returns a \code{list}, parallelizable.
# #' @param tree \code{TreeMan} object
# #' @param txnyms vector of taxonomic names
# #' @param parallel logical, make parallel? #' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
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