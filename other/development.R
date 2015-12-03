# TODO:
# -- viz function
# -- read newick
# -- get functions
# -- map and pin functions


library (treeman)
tree <- randTree (2)
viz (tree)


viz <- function (tree) {
  # y is a function of prenode's y and its own pd
  #  -- calculate the proportion of (pd + span) of the node
  #  -- cut the 2 space by the proportion
  #  -- - 1 from cut space
  #  -- the bigger the difference in y, the greater the y values
  get_pnts <- function (node, y, pnts) {
    pstnds <- node$postnode
    for (pstnd in pstnds) {
      pstnd <- tree@nodelist[[pstnd]]
      x <- pstnd$predist
      prp_y <- (2 * (pstnd$span + pstnd$pd)/node$pd) - 1
      y <- prp_y + y
      pnts <- rbind (pnts, data.frame (node=pstnd$id,
                                       x=x, y=y))
      pnts <- get_pnts (pstnd, y, pnts)
    }
    pnts
  }
  # start with root node
  # TODO: handle unrooted tree
  pnts <- data.frame (node=tree@root, x=0, y=0)
  root_node <- tree@nodelist[[tree@root]]
  pnts <- get_pnts (root_node, y=0, pnts=pnts)
  # add 10% to min y limit for node label
  y_lmts <- c (min(pnts$y) + (min(pnts$y)*.1), max(pnts$y))
  plot.default (x=pnts$x, y=pnts$y, col='black', pch=19, yaxt='n', ylab='',
                xlab='Branch span', bty='n', ylim=y_lmts)
  text (x=pnts$x, y=pnts$y, labels=pnts$node, pos=1)
  # draw lines
  for (i in 2:nrow (pnts)) {
    prenode <- tree@nodelist[[pnts$node[i]]]$prenode
    ind <- c (i, which (pnts$node == prenode))
    lines (x=pnts$x[ind], y=pnts$y[ind])
  }
}