randTree <- function (n) {
  .node <- function (n_left, id, span, prenode, predist, nodelist) {
    predist <- predist + span
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
    if (n_left >= 1) {
      n_left[2] <- ceiling (runif (min=0, max=n_left, n=1))
      n_left[1] <- abs (n_left[1] - n_left[2])
      for (n in n_left) {
        new_id <- paste0 ('n', length (nodelist) + 1)
        postnode <- c (postnode, new_id)
        new_span <- runif (min=0, max=1, n=1)
        pd <- new_span + pd
        new_prenode <- id
        new_predist <- predist + new_span
        nodelist <- .node (n, new_id, new_span, new_prenode, new_predist, nodelist)
        children <- c (children, nodelist[[new_id]]$children)
      }
      children <- c (postnode, children)
    }
    nodelist[[id]]$children <- children
    nodelist[[id]]$postnode <- postnode
    nodelist[[id]]$pd <- pd
    nodelist
  }
  # init empty nodelist
  nodelist <- list ()
  n_left <- (n * 2) - 1
  # generate root node
  id <- paste0 ('n', 1)
  predist <- span <- 0
  prenode <- c ()
  # gen nodelist
  nodelist <- .node (n_left, id, span, prenode, predist, nodelist)
  # init new tree object
  tree <- new ('NodeList', nodelist=nodelist, root='n1')
  .update (tree)
}