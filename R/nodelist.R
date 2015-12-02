#' @name NodeList
#' @title NodeList Class
#' @description \code{MoreTreeTool}'s S4 Class for representing phylogenetic trees.
#' A \code{NodeList} object holds a list of \code{Node} Reference Class Objects, each
#' of which points to connecting nodes. Changing one node changes all other nodes in the
#' tree. This allows for more efficient tree building, modelling and manipulation.
#' For example, adding a new node to a tree using a \code{phylo} class requires the
#' regeneration of the tree edge matrix. Whereas for a \code{NodeList} class the adding
#' of a new tip is computationally equivalent to adding a new element to a list.
#' @details Currently available methods:
#' \itemize{
#'   \item \code{tips()}, list all tips
#'   \item \code{nodes()}, list all internal nodes
#'   \item \code{nTips()}, count all tips
#'   \item \code{nNodes()}, count all internal nodes
#'   \item \code{rootNode()}, return root node ID
#'   \item \code{[[]]}, extract \code{Node}
#'   \item \code{pd()}, get total branch length
#'   \item \code{age()}, get max root to tip distance
#'   \item \code{ultrmtrc()}, is ultrametric T/F
#'   \item \code{plytms()}, is polytomous T/F
#'   \item \code{extant()}, return extant tips
#'   \item \code{extinct()}, return extinct tips
#'   \item \code{setTol()}, set tolerance (default 1e-8)
#' }
#' 
#' See examples for these methods in use.
#' @seealso
#' \code{\link{Node}}
#' @exportClass NodeList
#' @examples
#' library (MoreTreeTools)
#' # phylo
#' tree <- rtree (5)
#' # convert to NodeList
#' tree_nodelist <- as (tree_phylo, 'NodeList')
#' # Current available methods
#' print (tree_nodelist)  # get all 
#' tips (tree_nodelist)  # return all tips IDs
#' nodes (tree_nodelist)  # return all internal node IDs
#' nTips (tree_nodelist)  # count all tips
#' nNodes (tree_nodelist)  # count all internal nodes
#' rootNode(tree_nodelist)  # identify root node
#' tree_nodelist[['t1']]  # return t1 node object
#' pd (tree_nodelist)  # return phylogenetic diversity
#' age (tree_nodelist)  # return age of tree
#' ultrmtrc (tree_nodelist)  # is ultrametric?
#' plytms (tree_nodelist)  # is polytomous?
#' extant (tree_nodelist)  # return all extant tip IDs
#' extinct (tree_nodelist)  # return all extinct tip IDs
#' tree_nodelist <- setTol (tree_nodelist, 10)  # reset tolerance, default 1e-8
#' # now tol is higher more tips will be classed as extant
#' extant (tree)
# TODO: create validity check
setClass ('NodeList', representation=representation (
  nodelist='list',       # list of Node objects
  nodes='vector',        # vector of Node ids that are internal nodes
  tips='vector',         # vector of Node ids that are tips
  age='numeric',         # numeric of max root to tip distance
  pd='numeric',          # numeric of total branch length of tree
  extant='vector',       # vector of Node ids of all tips with 0 age
  extinct='vector',      # vector of Node ids of all tips with age > 0
  brnchlngth='logical',  # logical, do nodes have span
  ultrmtrc='logical',    # logical, do all tips end at 0
  plytms='logical',      # logical, is tree bifurcating
  tol='numeric',         # numeric of tolerance for determining extant
  root='character'),     # character of Node id of root
  prototype=prototype (tol=1e-8))

# Display methods
setMethod ('as.character', c('x'='NodeList'),
           function(x) {
             paste0 ('Tree (NodeList Object) of [', length (x@tips),'] tips')
           })
setMethod ('show', 'NodeList',
           function(object){
             print (object)
           })
setMethod ('str', c('object'='NodeList'),
           function (object, max.level=2L, ...) {
             if (is.na (max.level)) {
               stop ('max.level must be numeric')
             }
             str@default (object, max.level=max.level, ...)
           })
setGeneric ('print')
setMethod ('print', c('x'='NodeList'),
           function(x){
             msg <- 'Tree (NodeList Object):\n'
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

# Manip methods
setMethod ('[[', c ('NodeList', 'character', 'missing'),
           function(x, i, j, ...) {
             x@nodelist[[i]]
           })
setGeneric ("tips<-", signature=c("x"),
            function (x, value) {
              standardGeneric("tips<-")
            })
setReplaceMethod ("tips", "NodeList",
                  function (x, value) {
                    if (any (duplicated (value))) {
                      stop ('Tip names must be unique')
                    }
                    old_tips <- x@tips
                    n <- length (old_tips)
                    if (n != length (value)) {
                      stop ('Incorrect number of replacement tips')
                    }
                    mis <- match (old_tips, names (x@nodelist))
                    for (i in 1:n) {
                      x@nodelist[[old_tips[i]]]$id <- value[i]
                    }
                    names (x@nodelist)[mis] <- value
                    .update (x)
                  })
setGeneric ("nodes<-", signature=c("x"),
            function (x, value) {
              standardGeneric("nodes<-")
            })
setReplaceMethod ("nodes", "NodeList",
                  function (x, value) {
                    if (any (duplicated (value))) {
                      stop ('Node names must be unique')
                    }
                    old_nodes <- x@nodes
                    n <- length (old_nodes)
                    if (n != length (value)) {
                      stop ('Incorrect number of replacement nodes')
                    }
                    mis <- match (old_nodes, names (x@nodelist))
                    for (i in 1:n) {
                      x@nodelist[[old_nodes[i]]]$id <- value[i]
                    }
                    names (x@nodelist)[mis] <- value
                    .update (x)
                  })
setGeneric ('.update', signature=c('x'),
            function(x) {
              genericFunction ('.update')
            })
setMethod ('.update', 'NodeList',
           function (x) {
             with_pstndes <- sapply (x@nodelist,
                                     function (x) length (x$postnode) == 0)
             x@tips <- names (with_pstndes)[with_pstndes]
             x@nodes <- names (with_pstndes)[!with_pstndes]
             x@brnchlngth <- all (sapply (x@nodelist, function (x) length (x$span) > 0))
             if (x@brnchlngth) {
               if (length (x@root) > 0) {
                 x@age <- max (sapply (x@nodelist, function (x) x$predist))
                 extant_is <- unlist (sapply (x@tips, function (i) {
                   (x@age - x@nodelist[[i]]$predist) <= x@tol}))
                 x@extant <- names (extant_is)[extant_is]
                 x@extinct <- x@tips[!x@tips %in% x@extant]
                 x@ultrmtrc <- all (x@tips %in% extant (x))
               }
               x@pd <- x@nodelist[[x@root]]$pd
             }
             x@plytms <- any (sapply (x@nodelist, function (x) length (x$postnode) > 2))
             initialize (x)
           })

# Accessor method
setGeneric ('setTol', signature=c('x', 'n'),
            function(x, n) {
              genericFunction ('setTol')
            })
setMethod ('setTol', c ('NodeList', 'numeric'),
           function(x, n){
             x@tol <- n
             .update (x)
           })

# Info methods (user friendly)
setGeneric ('tips', signature=c('x'),
            function(x) {
              genericFunction ('tips')
            })
setMethod ('tips', 'NodeList',
           function(x){
             x@tips
           })
setGeneric ('nodes', signature=c('x'),
            function(x) {
              genericFunction ('nodes')
            })
setMethod ('nodes', 'NodeList',
           function(x){
             x@nodes
           })
setGeneric ('nTips', signature=c('x'),
            function(x) {
              genericFunction ('nTips')
            })
setMethod ('nTips', 'NodeList',
           function(x){
             length (x@tips)
           })
setGeneric ('plytms', signature=c('x'),
            function(x) {
              genericFunction ('plytms')
            }) 
setMethod ('plytms', 'NodeList',
           function(x) {
             x@plytms
           })
setGeneric ('ultrmtrc', signature=c('x'),
            function(x) {
              genericFunction ('ultrmtrc')
            }) 
setMethod ('ultrmtrc', 'NodeList',
           function(x) {
             x@ultrmtrc
           })
setGeneric ('extant', signature=c('x'),
            function(x) {
              genericFunction ('extant')
            }) 
setMethod ('extant', 'NodeList',
           function(x) {
             if (length (x@extant) == 0) {
               return (NULL)
             }
             x@extant
           })
setGeneric ('extinct', signature=c('x'),
            function(x) {
              genericFunction ('extinct')
            }) 
setMethod ('extinct', 'NodeList',
           function(x) {
             if (length (x@extinct) == 0) {
               return (NULL)
             }
             x@extinct
           })
setGeneric ('nNodes', signature=c('x'),
            function(x) {
              genericFunction ('nNodes')
            })
setMethod ('nNodes', 'NodeList',
           function(x){
             length (x@nodes)
           })
setGeneric ('rootNode', signature=c('x'),
            function(x) {
              genericFunction ('rootNode')
            })
setMethod ('rootNode', 'NodeList',
           function(x) {
             if (length (x@root) == 0) {
               return (NA)
             }
             x@root
           })
setGeneric ('age', signature=c('x'),
            function(x) {
              genericFunction ('age')
            })
setMethod ('age', 'NodeList',
           function(x){
             if (length (x@age) == 0) {
               return (NA)
             }
             x@age
           })
setGeneric ('pd', signature=c('x'),
            function(x) {
              genericFunction ('pd')
            })
setMethod ('pd', 'NodeList',
           function(x){
             if (length (x@pd) == 0) {
               return (NA)
             }
             x@pd
           })