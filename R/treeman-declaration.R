#' @name TreeMan
#' @title TreeMan Class
#' @description S4 Class for representing phylogenetic trees as a list of nodes.
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
#'    \item \code{taxonym}: name of taxonomic clade (optional)
#'    \item \code{span}: length of the preceding branch
#'    \item \code{pre}: ID of the preceding node
#'    \item \code{post}: IDs of the connecting nodes
#'    \item \code{children}: descending tip IDs
#'    \item \code{pd}: phylogenetic diversity represented by node
#'    \item \code{predist}: pre distance(distance to root if rooted or
#'    most distal tip if unrooted)
#' }
#' These data slots are updated whenever a node is modified, added or removed.
#' 
#' Currently available methods:
#' \itemize{
#'   \item \code{tips()}: list all tips
#'   \item \code{nodes()}: list all internal nodes
#'   \item \code{nTips()}: count all tips
#'   \item \code{nNodes()}: count all internal nodes
#'   \item \code{rootNode()}: return root node ID, NULL if unrooted
#'   \item \code{[[]]}: extract \code{Node}
#'   \item \code{pd()}: get total branch length of tree
#'   \item \code{age()}: get max root to tip distance
#'   \item \code{ultrmtrc()}: is ultrametric T/F
#'   \item \code{plytms()}: is polytomous T/F
#'   \item \code{extant()}: return extant tips
#'   \item \code{extinct()}: return extinct tips
#'   \item \code{setTol()}: set tolerance(default 1e-8)
#' }
#' 
#' See below in 'Examples' for these methods in use.
#' @seealso
#' \code{\link{randTree}}
#' @exportClass TreeMan
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
#' node['age']  # .... nchildren, pd, etc.

#TODO: modify this to allow user-defined Node slot
#TODO: check for missing children or no pd
essential_node_slots <- c('id')
valid_node_slots <- c('id', 'taxonym', 'span', 'prid',
                      'ptid', 'children', 'prdst', 'pd')

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
              stop(paste0(i, ' not in tree'))
            }
            .newNode(x, i)
          })
setMethod('[', c('TreeMan', 'character'),
          function(x, i) {
            if(!i %in% slotNames(x)) {
              stop(paste0(i, '  not in tree'))
            }
            slot(x, i)
          })