
#' @name pstMnp
#' @title Update prinds and tinds
#' @description Return tree with updated slots.
#' @details This function is automatically run. Only run, if you
#' are creating yor own functions to add and remove elements of the
#' \code{ndlst}.
#' @param tree \code{TreeMan} object
#' @seealso
#' \code{\link{updateSlts}}, \code{\link{addNdmtrx}},
#' \code{\link{getAge}}
#' @export
pstMnp <- function(tree) {
  # after any adding or removing of tips and nodes,
  # these slots MUST be updated to ensure full functionality
  tree@tinds <- .getTinds(tree@ndlst)
  tree@prinds <- .getPrinds(tree@ndlst)
  tree
}

#' @name updateSlts
#' @title Update tree slots after manipulation
#' @description Return tree with updated slots.
#' @details Tree slots in the \code{TreeMan} object are usually automatically updated.
#' For certain single node manipulations they are not. Run this
#' function to update the slots.
#' @param tree \code{TreeMan} object
#' @seealso
#' \code{\link{addNdmtrx}}, \code{\link{getAge}}
#' @export
updateSlts <- function(tree) {
  # Update the slots for a tree
  wo_pstndes <- vapply(tree@ndlst,
                       function(n) length(n[['ptid']]) == 0,
                       logical(1))
  tree@tips <- sort(names(wo_pstndes)[wo_pstndes])
  tree@ntips <- length(tree@tips)
  tree@nds <- sort(names(wo_pstndes)[!wo_pstndes])
  tree@nnds <- length(tree@nds)
  tree@all <- names(tree@ndlst)
  tree@nall <- length(tree@all)
  tree@wtxnyms <- any(vapply(tree@ndlst, function(n) !is.null(n[['txnym']]),
                             logical(1)))
  spns <- vapply(tree@ndlst, function(n) n[['spn']], numeric(1))
  tree@wspn <- any(spns > 0)
  if(tree@wspn) {
    tree@pd <- sum(vapply(tree@ndlst, function(n) n[['spn']], numeric(1)))
  } else {
    tree@pd <- numeric()
  }
  tree@ply <- any(vapply(tree@ndlst, function(n) length(n[['ptid']]) > 2,
                         logical(1)))
  tree@updtd <- TRUE
  initialize(tree)
}

#' @name addNdmtrx
#' @title Add node matrix to a tree
#' @description Return tree with node matrix added.
#' @details The node matrix makes 'enquiry'-type computations faster:
#' determining node ages, number of descendants etc. But it takes up
#' large amounts of memory and has no impact on adding or removing tips.
#' Note, trees with the node matrix can not be written to disk using the
#' 'serialization format' i.e. with \code{save} or \code{saveRDS}.
#' The matrix is generated with bigmemory's `as.big.matrix()`.
#' @param tree \code{TreeMan} object
#' @param shared T/F, should the bigmatrix be shared? See bigmemory documentation.
#' @param ... \code{as.big.matrix()} additional arguments
#' @seealso
#' \code{\link{updateSlts}}, \code{\link{rmNdmtrx}},
#' \url{https://cran.r-project.org/package=bigmemory}
#' @export
#' @examples
#' # library(treeman)
#' tree <- randTree(10, wndmtrx=FALSE)
#' summary(tree)
#' tree <- addNdmtrx(tree)
#' summary(tree)
addNdmtrx <- function(tree, shared=FALSE, ...) {
  if(tree@ntips < 3) {
    stop('Too small for node matrix.')
  }
  if(!checkNdlst(tree@ndlst, tree@root)) {
    stop('Invalid tree')
  }
  if(is.null(tree@ndmtrx)) {
    # generate ndmtrx
    tree@ndmtrx <- .getNdmtrxFrmLst(tree@ndlst, shared=shared, ...)
  }
  tree
}

#' @name rmNdmtrx
#' @title Remove node matrix
#' @description Return tree with memory heavy node matrix removed.
#' @details Potential uses: reduce memory load of a tree,
#' save tree using serialization methods.
#' @param tree \code{TreeMan} object
#' @seealso
#' \code{\link{addNdmtrx}}
#' @export
#' @examples
#' # library(treeman)
#' tree <- randTree(10)
#' summary(tree)
#' tree <- rmNdmtrx(tree)
#' summary(tree)
rmNdmtrx <- function(tree) {
  tree@ndmtrx <- NULL
  tree
}