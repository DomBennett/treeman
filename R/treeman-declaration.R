#TODO: modify this to allow user-defined Node slot
#TODO: check for missing kids or no pd
essential_node_slots <- c('id')
valid_node_slots <- c('id', 'txnym', 'span', 'prid',
                      'ptid', 'kids', 'prdst', 'pd')

.checkTreeMan <- function(object) {
  .check <- function(node) {
    if(!all(essential_node_slots %in% names(node))) {
      return(FALSE)
    }
    if(!all(names(node) %in% valid_node_slots)) {
      invlds <- names(node)[!names(node) %in% valid_node_slots]
      invlds <- paste(invlds, collapse="`, `")
      invlds <- paste0('[`', invlds, '`]')
      warning(paste0('Node [', node[['id']],
                     '] contains the following invalid node slots:\n',
                     invlds))
    }
    test_1 <- node[['id']] %in% nodes
    test_2 <- is.null(node[['prid']]) || (node[['prid']] %in% nodes)
    test_3 <- is.null(node[['prid']]) || all(node[['ptid']] %in% nodes)
    if(test_1 & test_2 & test_3) {
      return(TRUE)
    }
    FALSE
  }
  nodes <- names(object@nodelist)
  node_checks <- unlist(lapply(object@nodelist, .check))
  if(!all(node_checks)) {
    msg <- 'These nodes are invalid:\n'
    bad <- which(!node_checks)
    for(i in bad[-length(bad)]) {
      msg <- paste0(msg, nodes[i], ', ')
    }
    msg <- paste0(msg, nodes[bad[length(bad)]], '\n\n')
    msg <- paste0(msg, 'They may be pointing to non-existent nodes in tree, their ID may not be a named element in `@nodelist` or they may have missing essential node elements.')
    cat(msg)
    return(FALSE)
  }
  TRUE
}

#' @name TreeMan
#' @title S4 class for representing phylogenetic trees as a list of nodes.
#' @slot nodelist list of nodes
#' @slot nodes vector of node ids that are internal nodes
#' @slot nnodes numeric of number of internal nodes in tree
#' @slot tips vector of node ids that are tips
#' @slot ntips numeric of number of internal nodes in tree
#' @slot all vector of all node ids
#' @slot nall numeric of number of all nodes in tree
#' @slot age numeric of max root to tip distance
#' @slot pd numeric of total branch length of tree
#' @slot ext vector of Node ids of all tips with 0 age
#' @slot exc vector of Node ids of all tips with age > 0
#' @slot wspn logical, do nodes have spans
#' @slot ultr logical, do all tips end at 0
#' @slot ply logical, is tree bifurcating
#' @slot tol numeric of tolerance for determining extant
#' @slot root character of Node id of root, if no root then empty character
#' @details
#' A \code{TreeMan} object holds a list of nodes. The idea of the \code{TreeMan}
#' class is to make adding and removing nodes as similar as possible to adding
#' and removing elements in a list. Note that internal nodes and tips are
#' both considered nodes. Trees can be unrooted and polytomous.
#' 
#' 
#' Each node within the \code{TreeMan} \code{nodelist} contains the following data slots:
#' \itemize{
#'    \item \code{id}: character string for the node ID
#'    \item \code{txnym}: name of taxonomic clade (optional)
#'    \item \code{span}: length of the preceding branch
#'    \item \code{prid}: IDs of the preceding nodes to the root
#'    \item \code{ptid}: IDs of the immediately connecting nodes
#'    \item \code{kids}: descending tip IDs
#'    \item \code{pd}: phylogenetic diversity represented by node
#'    \item \code{prdst}: pre distance(distance to root if rooted or
#'    most distal tip if unrooted)
#' }
#' These data slots are updated whenever a node is modified, added or removed.
#' 
#' See below in 'Examples' for these methods in use.
#' @seealso
#' \code{\link{randTree}}, \code{\link{Node}}
#' @examples
#' library(treeman)
#' # Generate random tree
#' tree <- randTree(10)
#' # Print to get basic stats
#' print(tree)
#' # Currently available methods
#' tree['tips']  # return all tips IDs
#' tree['nodes']  # return all internal node IDs
#' tree['ntips']  # count all tips
#' tree['nnodes']  # count all internal nodes
#' tree['root']  # identify root node
#' tree[['t1']]  # return t1 node object
#' tree['pd']  # return phylogenetic diversity
#' tree['age']  # return age of tree
#' tree['ultr']  # is ultrametric?
#' tree['ply']  # is polytomous?
#' tree['ext']  # return all extant tip IDs
#' tree['exc']  # return all extinct tip IDs
#' tree <- setTol(tree, 10)  # reset tolerance, default 1e-8
#' # now tol is higher more tips will be classed as extant
#' tree['ext']
#' # Because all nodes are lists with metadata we can readily
#' #  get specific information on nodes of interest
#' node <- tree[['n2']]
#' print(node)
#' # And then use the same syntax for the tree
#' node['age']  # .... nkids, pd, etc.
#' @exportClass TreeMan
setClass('TreeMan', representation=representation(
  nodelist='list',       # list of Node objects
  nodes='vector',        # vector of Node ids that are internal nodes
  nnodes='numeric',      # numeric of number of internal nodes in tree
  tips='vector',         # vector of Node ids that are tips
  ntips='numeric',       # numeric of number of internal nodes in tree
  all='vector',          # vector of all Node ids
  nall='numeric',        # numeric of number of all nodes in tree
  age='numeric',         # numeric of max root to tip distance
  pd='numeric',          # numeric of total branch length of tree
  ext='vector',          # vector of Node ids of all tips with 0 age
  exc='vector',          # vector of Node ids of all tips with age > 0
  wspn='logical',        # logical, do nodes have spans
  ultr='logical',        # logical, do all tips end at 0
  ply='logical',         # logical, is tree bifurcating
  tol='numeric',         # numeric of tolerance for determining extant
  root='character'),     # character of Node id of root, if no root then empty character
  prototype=prototype(tol=1e-8), validity=.checkTreeMan)

# Accessor methods
setMethod('[[', c('TreeMan', 'character'),
          function(x, i) {
            if(!i %in% names(x@nodelist)) {
              srch_trm <- gsub(' ', '_', i)  # usual mistake
              pssbls <- which(agrepl(srch_trm, names(x@nodelist), ignore.case=TRUE,
                                     max.distance=0.25))
              pssbls <- names(x@nodelist)[pssbls]
              if(length(pssbls) > 0 & length(pssbls) < 50) {
                msg <- paste0("Can't find [", i, "]. Did you mean ....\n")
                for(p in pssbls) {
                  msg <- paste0(msg, '"', p, '"\n')
                }
                msg <- paste0(msg, "?\n")
              } else {
                msg <- paste0("Can't find [", i, "] in tree.")
              }
              stop(msg)
            }
            .newNode(x, i)
          })
setMethod('[', c('TreeMan', 'character'),
          function(x, i) {
            slt_nms <- slotNames(x)
            slt_nms <- slt_nms[slt_nms != 'ndlst']
            if(!i %in% slt_nms) {
              slt_nms <- paste0(slt_nms, collapse=', ')
              stop(paste0('`', i, '` not a tree slot. Available slots: ', slt_nms))
            }
            slot(x, i)
          })

# display methods
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

#' @name viz
#' @title Visualise a tree
#' @param tree \code{TreeMan}
#' @param taxonyms T/F, show the nodes with their taxonyms instead of IDs
#' @description Crude plot of \code{TreeMan} objects
#' @export
setGeneric("viz", signature=c("tree", "taxonyms"),
           function(tree, taxonyms=FALSE) {
             standardGeneric("viz")
           })
#' @exportMethod viz
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
                tree@nodelist[[i]][['pd']] <- length(tree@nodelist[[i]][['kids']])
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
              prenode <- tree@nodelist[[pnts$node[i]]][['prid']][1]
              ind <- c(i, which(pnts$node == prenode))
              lines(x=pnts$x[ind], y=pnts$y[ind])
            }
          })