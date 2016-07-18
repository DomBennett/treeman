#' @name randTree
#' @title Generate a random tree
#' @description Returns a random \code{TreeMan} tree with \code{n}
#' tips.
#' @details Equivalent to \code{ape}'s \code{rtree()} but returns a
#' \code{TreeMan} tree. Tree is always rooted and bifurcating.
#' @param n number of tips, integer, must be 3 or greater
#' @seealso
#' \code{\link{TreeMan-class}}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(5)
randTree <- function(n) {
  # Return a random tree based on a broken-stick model
  .add <- function(i) {
    nd <- vector("list", length=4)
    names(nd) <- c('id', 'ptid', 'prid', 'spn')
    nd[['id']] <- ids[i]
    nd[['spn']] <- spns[i]
    nd[['prid']] <- ids[prids[i]]
    ptids <- ids[prids == i]
    nd[['ptid']] <- ptids[!is.na(ptids)]
    nd
  }
  if(n < 3) {
    stop("n is too small")
  }
  nnds <- n + (n - 1)
  prids <- rep(NA, nnds)
  tmp_ids <- 1:nnds
  pool <- rep(2, nnds)
  pool[1] <- 0
  prids[1:3] <- 1
  for(i in seq(4, nnds, 2)) {
    j <- sample(1:(i-1), 1, prob=pool[1:(i-1)])
    pool[j] <- pool[j] - 1
    prids[i] <- prids[i + 1] <- tmp_ids[j]
  }
  spns <- c(0, runif(nnds-1, 0, 1))
  tids <- paste0('t', 1:n)
  nids <- paste0('n', 1:(n-1))
  ids <- rep(NA, nnds)
  ids[tmp_ids %in% prids] <- nids
  sum(!tmp_ids %in% prids)
  ids[!tmp_ids %in% prids] <- tids
  ndlst <- lapply(1:nnds, .add)
  names(ndlst) <- ids
  # init new tree object
  tree <- new('TreeMan', ndlst=ndlst, root='n1')
  updateTree(tree)
}

blncdTree <- function(...) {
  cat('This function is in progress.... \n')
}

imblncdTree <- function(...) {
  cat('This function is in progress.... \n')
}