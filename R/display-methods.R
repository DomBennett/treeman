# TreeMen display methods
setMethod('as.character', c('x'='TreeMen'),
          function(x) {
            paste0('TreeMen Object of [', x@ntrees,'] trees')
          })
setMethod('show', 'TreeMen',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
setMethod('str', c('object'='TreeMen'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
setMethod('print', c('x'='TreeMen'),
          function(x){
            msg <- 'Trees (TreeMen Object):\n'
            msg <- paste0(msg, '  + ', x@ntrees, ' trees\n')
            msg <- paste0(msg, '  + ', x@ntips, ' tips\n')
            cat(msg)
          })


# TreeMan display methods
setMethod('as.character', c('x'='TreeMan'),
           function(x) {
             paste0('TreeMan Object of [', length(x@tips),'] tips')
           })
setMethod('show', 'TreeMan',
           function(object){
             msg <- as.character(object)
             cat(msg)
           })
setMethod('str', c('object'='TreeMan'),
           function(object, max.level=2L, ...) {
             if(is.na(max.level)) {
               stop('max.level must be numeric')
             }
             str@default(object, max.level=max.level, ...)
           })
setMethod('print', c('x'='TreeMan'),
           function(x){
             msg <- 'Tree (TreeMan Object):\n'
             msg <- paste0(msg, '  + ', x@ntips, ' tips\n')
             msg <- paste0(msg, '  + ', x@nnodes, ' internal nodes\n')
             if(x@ply) {
               msg <- paste0(msg, '  + Polytomous\n')
             } else {
               msg <- paste0(msg, '  + Binary\n')
             }
             if(is.na(x@root)) {
               if(x@wspn) {
                 msg <- paste0(msg, '  + Unrooted and without node spans\n')
               } else {
                 msg <- paste0(msg, '  + Unrooted, with node spans\n')
                 msg <- paste0(msg, '  + PD ', signif(x@pd, 3), '\n')
               }
             } else {
               if(x@wspn) {
                 msg <- paste0(msg, '  + Age ', signif(x@age, 3), '\n')
                 msg <- paste0(msg, '  + PD ', signif(x@pd, 3), '\n')
                 if(x@ultr) {
                   msg <- paste0(msg, '  + Ultrametric (all tips are extant)\n')
                 } else {
                   msg <- paste0(msg, '  + Not ultrametric (with extinct tips)\n')
                 }
               } else {
                 msg <- paste0(msg, '  + Without node spans\n')
               }
               msg <- paste0(msg, '  + Root node is \"', x@root, '\"\n')
             }
             cat(msg)
           })
setGeneric("viz", signature=c("tree", "taxonyms"),
            function(tree, taxonyms=FALSE) {
              standardGeneric("viz")
            })
setMethod('viz', 'TreeMan',
           function(tree, taxonyms){
             get_pnts <- function(node, y, pnts) {
               pstids <- node[['ptid']]
               low_y_diff <- -node[['pd']]/2
               high_y_diff <- node[['pd']]/2
               y_diffs <- seq(from=low_y_diff, to=high_y_diff,
                              length.out=length(pstids))
               counter <- 1
               for(pstid in pstids) {
                 pstnd <- tree@nodelist[[pstid]]
                 pstnd_x <- pstnd[['prdst']]
                 pstnd_y <- y + y_diffs[counter]
                 pnts <- rbind(pnts, data.frame(node=pstid,
                                                  x=pstnd_x, y=pstnd_y))
                 pnts <- get_pnts(pstnd, pstnd_y, pnts)
                 counter <- counter + 1
               }
               pnts
             }
             if(!tree@wspn) {
               # TODO: switch to setNodesSpan
               for(i in 1:length(tree@nodelist)) {
                 tree@nodelist[[i]][['span']] <- 1
                 tree@nodelist[[i]][['pd']] <- length(tree@nodelist[[i]][['children']])
                 prids <- getNodePrid(tree, tree@nodelist[[i]][['id']])
                 tree@nodelist[[i]][['prdst']] <- length(prids)
               }
               tree@pd <- length(tree@nodelist) - 1
             }
             # start with root node
             # TODO: handle unrooted tree
             pnts <- data.frame(node=tree@root, x=0, y=tree@pd, stringsAsFactors=FALSE)
             root_node <- tree@nodelist[[tree@root]]
             pnts <- get_pnts(root_node, y=tree@pd, pnts=pnts)
             # add 10% to min y limit for node label
             min_y <- abs(min(pnts$y))
             min_y <- min_y + (min_y*.1)
             min_y <- ifelse(min(pnts$y) > 0, min_y, -1*min_y)
             y_lmts <- c(min_y, max(pnts$y))
             plot.default(x=pnts$x, y=pnts$y, col='black', pch=19, yaxt='n', ylab='',
                           xlab='', bty='n', ylim=y_lmts)
             if(taxonyms) {
               text(x=pnts$x, y=pnts$y,
                     labels=sapply(pnts$node, function(n) tree@nodelist[[n]][['taxonym']]),
                     pos=1)
             } else {
               text(x=pnts$x, y=pnts$y, labels=pnts$node, pos=1)
             }
             # draw lines
             for(i in 2:nrow (pnts)) {
               prenode <- tree@nodelist[[pnts$node[i]]][['prid']]
               ind <- c(i, which(pnts$node == prenode))
               lines(x=pnts$x[ind], y=pnts$y[ind])
             }
           })