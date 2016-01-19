#' @name randTree
#' @title Generate a random tree
#' @description Returns a random \code{TreeMan} tree with \code{n}
#' tips.
#' @details Equivalent to \code{ape}'s \code{rtree()} but returns a
#' \code{TreeMan} tree. Tree is always rooted and bifurcating.
#' @param n number of tips, integer, must be greater than 1
#' @seealso
#' \code{\link{TreeMan}}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(5)

randTree <- function (n) {
  .node <- function (n_left, id, span, prenode, predist, nodelist) {
    postnode <- children <- c ()
    pd <- 0
    node <- list ('id'=id,
                  'span'=span,
                  'prenode'=prenode,
                  'postnode'=postnode,
                  'children'=children,
                  'pd'=pd,
                  'predist'=predist)
    nodelist[[id]] <- node
    # if there are enough ns left to have children
    n_left <- n_left - 1
    if (n_left > 1) {
      nls <- n_left
      if (n_left == 2) {
        nls <- c (1, 1)
      } else {
        # split must be binary
        splits <- seq (from=1, to=((n_left)-1)/2, by=2)
        nls[2] <- sample (splits, 1)
        nls[1] <- n_left - nls[2]
      }
      for (nl in nls) {
        if (nl == 1) {
          ntips <<- ntips + 1
          new_id <- paste0 ('t', ntips)
          children <- c (new_id, children)
        } else {
          nnodes <<- nnodes + 1
          new_id <- paste0 ('n', nnodes)
        }
        postnode <- c (postnode, new_id)
        new_span <- runif (min=0, max=1, n=1)
        new_prenode <- id
        new_predist <- predist + new_span
        nodelist <- .node (nl, new_id, new_span,
                           new_prenode, new_predist, nodelist)
        children <- c (children, nodelist[[new_id]]$children)
        pd <- pd + new_span + nodelist[[new_id]]$pd
      }
    }
    nodelist[[id]]$children <- children
    nodelist[[id]]$postnode <- postnode
    nodelist[[id]]$pd <- pd
    nodelist
  }
  if (n < 2) {
    stop('`n` must be greater than 1')
  }
  # init empty nodelist
  ntips <- 0
  nnodes <- 1
  nodelist <- list ()
  n_left <- (n - 1) + n
  # generate root node
  id <- paste0 ('n', 1)
  predist <- span <- 0
  prenode <- c ()
  # gen nodelist
  nodelist <- .node (n_left, id, span, prenode, predist, nodelist)
  # init new tree object
  tree <- new ('TreeMan', nodelist=nodelist, root='n1')
  .update (tree)
}

#TODO balancedTree