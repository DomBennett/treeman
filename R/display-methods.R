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
setGeneric ("viz", signature=c("tree", "taxonyms"),
            function (tree, taxonyms=FALSE) {
              standardGeneric("viz")
            })
setMethod ('viz', 'TreeMan',
           function(tree, taxonyms){
             get_pnts <- function (node, y, pnts) {
               pstnds <- node$postnode
               low_y_diff <- -node$pd/2
               high_y_diff <- node$pd/2
               y_diffs <- seq(from=low_y_diff, to=high_y_diff,
                              length.out=length(pstnds))
               counter <- 1
               for (pstnd in pstnds) {
                 pstnd <- tree@nodelist[[pstnd]]
                 pstnd_x <- pstnd$predist
                 pstnd_y <- y + y_diffs[counter]
                 pnts <- rbind (pnts, data.frame (node=pstnd$id,
                                                  x=pstnd_x, y=pstnd_y))
                 pnts <- get_pnts (pstnd, pstnd_y, pnts)
                 counter <- counter + 1
               }
               pnts
             }
             # start with root node
             # TODO: handle unrooted tree
             pnts <- data.frame (node=tree@root, x=0, y=tree@pd, stringsAsFactors=FALSE)
             root_node <- tree@nodelist[[tree@root]]
             pnts <- get_pnts (root_node, y=tree@pd, pnts=pnts)
             # add 20% to min y limit for node label
             y_lmts <- c (min(pnts$y) - (min(pnts$y)*.2), max(pnts$y))
             plot.default (x=pnts$x, y=pnts$y, col='black', pch=19, yaxt='n', ylab='',
                           xlab='', bty='n', ylim=y_lmts)
             if(taxonyms) {
               text (x=pnts$x, y=pnts$y,
                     labels=sapply(pnts$node, function(n) tree@nodelist[[n]]$taxonym),
                     pos=1)
             } else {
               text (x=pnts$x, y=pnts$y, labels=pnts$node, pos=1)
             }
             # draw lines
             for (i in 2:nrow (pnts)) {
               prenode <- tree@nodelist[[pnts$node[i]]]$prenode
               ind <- c (i, which (pnts$node == prenode))
               lines (x=pnts$x[ind], y=pnts$y[ind])
             }
           })