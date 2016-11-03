#' @name randTree
#' @title Generate a random tree
#' @description Returns a random \code{TreeMan} tree with \code{n}
#' tips.
#' @details Equivalent to \code{ape}'s \code{rtree()} but returns a
#' \code{TreeMan} tree. Tree is always rooted and bifurcating.
#' @param n number of tips, integer, must be 3 or greater
#' @param wndmtrx T/F add node matrix? Default TRUE.
#' @param parallel T/F run in parallel? Default FALSE.
#' @seealso
#' \code{\link{TreeMan-class}}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(5)
randTree <- function(n, wndmtrx=TRUE, parallel=FALSE) {
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
  tree <- new('TreeMan', ndlst=ndlst, root='n1', wtxnyms=FALSE,
              ndmtrx=NULL, prinds=prinds, tinds=tinds)
  tree <- updateSlts(tree)
  if(wndmtrx) {
    tree <- addNdmtrx(tree)
  }
  tree
}

blncdTree <- function(...) {
  cat('This function is in progress.... \n')
}

imblncdTree <- function(...) {
  cat('This function is in progress.... \n')
}