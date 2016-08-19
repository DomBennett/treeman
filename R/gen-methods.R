#' @name randTree
#' @title Generate a random tree
#' @description Returns a random \code{TreeMan} tree with \code{n}
#' tips.
#' @details Equivalent to \code{ape}'s \code{rtree()} but returns a
#' \code{TreeMan} tree. Tree is always rooted and bifurcating.
#' @param n number of tips, integer, must be 3 or greater
#' @param update T/F update tree slots after generation? Default TRUE.
#' @param parallel T/F run in parallel? Default FALSE.
#' @seealso
#' \code{\link{TreeMan-class}}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(5)
randTree <- function(n, update=TRUE, parallel=FALSE) {
  # Return a random tree based on a broken-stick model
  .add <- function(i) {
    nd <- vector("list", length=4)
    names(nd) <- c('id', 'ptid', 'prid', 'spn')
    nd[['id']] <- ids[i]
    nd[['spn']] <- spns[i]
    nd[['prid']] <- ids[prinds[i]]
    nd[['ptid']] <- ptids[ptnds_pool == i]
    nd
  }
  if(n < 3) {
    stop("`n` is too small")
  }
  nnds <- n + (n - 1)
  prinds <- rep(NA, nnds)
  intrnls <- seq(2, (nnds-2), 2)
  # randomise intrnls
  intrnls <- intrnls +
    sample(0:1, size=length(intrnls), replace=TRUE)
  # create prids vector
  prinds[1:3] <- 1
  prinds[4:nnds] <- rep(intrnls, each=2)
  # random numbers for spans
  spns <- c(0, runif(nnds-1, 0, 1))
  ids <- rep(NA, nnds)
  tinds <- which(!1:nnds %in% prinds)
  ids[tinds] <- paste0('t', 1:n)
  ids[1:nnds %in% prinds] <- paste0('n', 1:(n-1))
  ptnds_pool <- prinds[-1]
  ptids <- ids[-1]
  ndlst <- plyr::mlply(.data=1:nnds, .fun=.add, .parallel=parallel)
  attr(ndlst, "split_labels") <- 
    attr(ndlst, "split_type") <- NULL
  names(ndlst) <- ids
  # init new tree object
  tree <- new('TreeMan', ndlst=ndlst, root='n1',
              ndmtrx=bigmemory::big.matrix(1,1),
              prinds=prinds, tinds=tinds)
  if(update) {
    tree <- updateTree(tree)
  } else {
    # init basic slots
    tree@updtd <- FALSE
    tree@tips <- paste0('t', 1:n)
    tree@ntips <- n
    tree@nds <- paste0('n', 1:(n-1))
    tree@nnds <- n - 1
    tree@all <- names(tree@ndlst)
    tree@nall <- nnds
    tree@wspn <- TRUE
  }
  tree
}

blncdTree <- function(...) {
  cat('This function is in progress.... \n')
}

imblncdTree <- function(...) {
  cat('This function is in progress.... \n')
}