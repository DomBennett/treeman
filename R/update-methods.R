#' @name updateTree
#' @title Update tree slots after manipulation
#' @description Return tree with updated node slots.
#' @details Tree slots in the \code{TreeMan} object are not automatically updated
#' to make mainipulations computationally faster. To update these slots, run this function.
#' @param tree \code{TreeMan} object
#' @seealso
#' \code{\link{downdateTree}}, \code{\link{getTreeAge}}
#' @export
updateTree <- function(tree) {
  if(!checkTreeMan(tree)) {
    stop('Invalid tree')
  }
  if(length(tree@tinds) == 0) {
    tree@tinds <- .getTinds(tree@ndlst)
  }
  if(length(tree@prinds) == 0) {
    tree@prinds <- .getPrinds(tree@ndlst)
  }
  if(is.null(tree@ndmtrx)) {
    # generate ndmtrx
    tree@ndmtrx <- .getNdmtrxFrmLst(tree@ndlst)
  }
  # Update the slots for a tree
  wo_pstndes <- sapply(tree@ndlst,
                       function(n) length(n[['ptid']]) == 0)
  tree@tips <- sort(names(wo_pstndes)[wo_pstndes])
  tree@ntips <- length(tree@tips)
  tree@nds <- sort(names(wo_pstndes)[!wo_pstndes])
  tree@nnds <- length(tree@nds)
  tree@all <- names(tree@ndlst)
  tree@nall <- length(tree@all)
  tree@wtxnyms <- any(sapply(tree@ndlst, function(n) !is.null(n[['txnym']])))
  spns <- sapply(tree@ndlst, function(n) n[['spn']])
  tree@wspn <- any(spns > 0)
  if(tree@wspn) {
    if(length(tree@root) > 0) {
      tip_prdsts <- .getNdsPrdstsFrmMtrx(tree@ndmtrx, tree@all,
                                         tree@tips, spns,
                                         parallel=FALSE,
                                         progress="none")
      tree@age <- max(tip_prdsts)
      extant_is <- (tree@age - tip_prdsts) <= tree@tol
      tree@ext <- names(extant_is)[extant_is]
      tree@exc <- tree@tips[!tree@tips %in% tree@ext]
      tree@ultr <- all(tree@tips %in% tree@ext)
    } else {
      tree@ext <- tree@exc <- vector()
      tree@ultr <- FALSE
      tree@age <- numeric()
    }
    tree@pd <- sum(sapply(tree@ndlst, function(n) n[['spn']]))
  } else {
    tree@age <- tree@pd <- numeric()
    tree@ext <- tree@ext <- vector()
    tree@ultr <- logical()
  }
  tree@ply <- any(sapply(tree@ndlst, function(n) length(n[['ptid']]) > 2))
  tree@updtd <- TRUE
  initialize(tree)
}

#' @name downdateTree
#' @title Downdate tree slots
#' @description Return tree with memory heavy node slots removed.
#' @details Potential uses: reduce memory load by downdating a tree,
#' force get-methods to use ndlst rather than ndmtrx.
#' @param tree \code{TreeMan} object
#' @seealso
#' \code{\link{updateTree}}
#' @export
#' @examples
#' # library(treeman)
#' tree <- randTree(10)
#' tree <- downdateTree(tree)
#' # summary(tree)  # running this will lead to an error
downdateTree <- function(tree) {
  tree@ndmtrx <- NULL
  tree@updtd <- FALSE
  tree@tinds <- vector("integer", length=0)
  tree@prinds <- vector("integer", length=0)
  tree
}