# TreeMan display methods

setMethod ('as.character', c('x'='TreeMan'),
           function(x) {
             paste0 ('TreeMan Object of [', length (x@tips),'] tips')
           })
setMethod ('show', 'TreeMan',
           function(object){
             print (object)
           })
setMethod ('str', c('object'='TreeMan'),
           function (object, max.level=2L, ...) {
             if (is.na (max.level)) {
               stop ('max.level must be numeric')
             }
             str@default (object, max.level=max.level, ...)
           })
setGeneric ('print')
setMethod ('print', c('x'='TreeMan'),
           function(x){
             msg <- 'Tree (TreeMan Object):\n'
             msg <- paste0 (msg, '  -- [', length (x@tips), '] tips\n')
             msg <- paste0 (msg, '  -- [', length (x@nodes), '] internal nodes\n')
             if (x@plytms) {
               msg <- paste0 (msg, '  -- Polytomous\n')
             } else {
               msg <- paste0 (msg, '  -- Binary\n')
             }
             if (is.na (rootNode (x))) {
               if (is.na (age (x))) {
                 msg <- paste0 (msg, '  -- Unrooted and without branch lengths\n')
               } else {
                 msg <- paste0 (msg, '  -- Unrooted, with branch lengths\n')
                 msg <- paste0 (msg, '  -- PD [', signif (x@pd, 3), ']\n')
               }
             } else {
               msg <- paste0 (msg, '  -- Age [', signif (x@age, 3), ']\n')
               msg <- paste0 (msg, '  -- PD [', signif (x@pd, 3), ']\n')
               if (x@ultrmtrc) {
                 msg <- paste0 (msg, '  -- Ultrametric (all tips are extant)\n')
               } else {
                 msg <- paste0 (msg, '  -- Not ultrametric (with extinct tips)\n')
               }
               msg <- paste0 (msg, '  -- Root node is ["', x@root, '"]\n')
             }
             cat (msg)
           })
setGeneric ("viz", signature=c("tree"),
            function (tree) {
              standardGeneric("viz")
            })
setMethod ('viz', 'TreeMan',
           function(tree){
             .viz (tree)
           })

# Testable function
.viz <- function (tree) {
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
                xlab='', bty='n', ylim=y_lmts)
  text (x=pnts$x, y=pnts$y, labels=pnts$node, pos=1)
  # draw lines
  for (i in 2:nrow (pnts)) {
    prenode <- tree@nodelist[[pnts$node[i]]]$prenode
    ind <- c (i, which (pnts$node == prenode))
    lines (x=pnts$x[ind], y=pnts$y[ind])
  }
}