#' @name updateTree
#' @title Update tree slots after manipulation
#' @description Return tree with updated node slots.
#' @details Tree slots in the \code{TreeMan} object are not automatically updated
#' to make mainipulations computationally faster. To update these slots, run this function.
#' @param tree \code{TreeMan} object
#' @seealso
#' \code{\link{getTreeAge}}
#' @export
#' @examples
#' # library(treeman)
#' # tree <- randTree(10)
#' # rmTip
#' # summary(tree)  # old information
#' # tree <- updateTree(tree)
#' # summary(tree)  # new information
updateTree <- function(tree) {
  if(!checkTreeMan(tree)) {
    stop('Invalid tree')
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
  tree@wspn <- any(sapply(tree@ndlst, function(n) n[['spn']] != 0))
  if(tree@wspn) {
    if(length(tree@root) > 0) {
      tip_prdsts <- getNdsPrdst(tree, tree@tips)
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
  initialize(tree)
}