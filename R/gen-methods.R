#' @name randTree
#' @title Generate a random tree
#' @description Returns a random \code{TreeMan} tree with \code{n}
#' tips.
#' @details Equivalent to \code{ape}'s \code{rtree()} but returns a
#' \code{TreeMan} tree. Tree is always rooted and bifurcating.
#' @param n number of tips, integer, must be 3 or greater
#' @param wndmtrx T/F add node matrix? Default FALSE.
#' @param parallel T/F run in parallel? Default FALSE.
#' @seealso
#' \code{\link{TreeMan-class}}, \code{\link{blncdTree}},
#' \code{\link{unblncdTree}}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(5)
randTree <- function(n, wndmtrx=FALSE, parallel=FALSE) {
  # Return a random tree based on a broken-stick model
  .randomPrinds <- function(n) {
    pool <- rep((1:(n-1)), each=2)
    res <- rep(NA, length(pool)+1)
    res[1] <- 1
    for(i in 2:length(res)) {
      pssbls <- which(i > pool)
      if(length(pssbls) == 1) {
        i_pool <- pssbls
      } else {
        i_pool <- sample(pssbls, 1)
      }
      res[i] <- pool[i_pool]
      pool[i_pool] <- NA
    }
    res
  }
  if(n < 3) {
    stop("`n` is too small")
  }
  prinds <- .randomPrinds(n)
  .cnstrctTree(n, prinds, wndmtrx=wndmtrx,
               parallel=parallel)
}

#' @name blncdTree
#' @title Generate a balanced tree
#' @description Returns a balanced \code{TreeMan} tree with \code{n}
#' tips.
#' @details Equivalent to \code{ape}'s \code{stree(type='balanced')} but returns a
#' \code{TreeMan} tree. Tree is always rooted and bifurcating.
#' @param n number of tips, integer, must be 3 or greater
#' @param wndmtrx T/F add node matrix? Default FALSE.
#' @param parallel T/F run in parallel? Default FALSE.
#' @seealso
#' \code{\link{TreeMan-class}}, \code{\link{randTree}},
#' \code{\link{unblncdTree}}
#' @export
#' @examples
#' library(treeman)
#' tree <- blncdTree(5)
blncdTree <- function(n, wndmtrx=FALSE, parallel=FALSE) {
  if(n < 3) {
    stop("`n` is too small")
  }
  prinds <- c(1, rep((1:(n-1)), each=2))
  .cnstrctTree(n, prinds, wndmtrx=wndmtrx,
               parallel=parallel)
}

#' @name unblncdTree
#' @title Generate an unbalanced tree
#' @description Returns an unbalanced \code{TreeMan} tree with \code{n}
#' tips.
#' @details Equivalent to \code{ape}'s \code{stree(type='left')} but returns a
#' \code{TreeMan} tree. Tree is always rooted and bifurcating.
#' @param n number of tips, integer, must be 3 or greater
#' @param wndmtrx T/F add node matrix? Default FALSE.
#' @param parallel T/F run in parallel? Default FALSE.
#' @seealso
#' \code{\link{TreeMan-class}}, \code{\link{randTree}},
#' \code{\link{blncdTree}}
#' @export
#' @examples
#' library(treeman)
#' tree <- unblncdTree(5)
unblncdTree <- function(n, wndmtrx=FALSE, parallel=FALSE) {
  if(n < 3) {
    stop("`n` is too small")
  }
  prinds <- c(1, 1:(n-1), 1:(n-1))
  .cnstrctTree(n, prinds, wndmtrx=wndmtrx,
               parallel=parallel)
}

.cnstrctTree <- function(n, prinds, wndmtrx, parallel) {
  .add <- function(i) {
    nd <- vector("list", length=4)
    names(nd) <- c('id', 'ptid', 'prid', 'spn')
    nd[['id']] <- ids[i]
    nd[['spn']] <- spns[i]
    nd[['prid']] <- ids[prinds[i]]
    nd[['ptid']] <- ptids[ptnds_pool == i]
    nd
  }
  nnds <- length(prinds)
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
  tree <- new('TreeMan', ndlst=ndlst, root='n1', wtxnyms=FALSE,
              ndmtrx=NULL, prinds=prinds, tinds=tinds)
  tree <- updateSlts(tree)
  if(wndmtrx) {
    tree <- addNdmtrx(tree)
  }
  tree
}

#' @name twoer
#' @title Generate a tree of two tips
#' @description Returns a \code{TreeMan} tree with two tips and a root.
#' @details Useful for building larger trees with \code{addClade()}.
#' Note, a node matrix cannot be added to a tree of two tips.
#' @param tids tip IDs
#' @param spns tip spans
#' @param rid root ID
#' @param root_spn root span
#' @seealso
#' \code{\link{TreeMan-class}}, \code{\link{randTree}}
#' @export
#' @examples
#' library(treeman)
#' tree <- twoer()
twoer <- function(tids=c('t1', 't2'), spns=c(1,1),
                  rid='root', root_spn=0) {
  ndlst <- list()
  ndlst[[rid]] <- list('id'=rid, 'prid'=rid,
                       'ptid'=tids[1:2], 'spn'=root_spn)
  ndlst[[tids[1]]] <- list('id'=tids[[1]], 'prid'=rid,
                           'ptid'=NULL, 'spn'=spns[1])
  ndlst[[tids[2]]] <- list('id'=tids[[2]], 'prid'=rid,
                           'ptid'=NULL, 'spn'=spns[2])
  prinds <- c(1, 1, 1)
  tinds <- c(2, 3)
  tree <- new('TreeMan', ndlst=ndlst, root='root', wtxnyms=FALSE,
              ndmtrx=NULL, prinds=prinds, tinds=tinds)
  updateSlts(tree)
}
