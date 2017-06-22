# roxygen imports
#' @import methods
#' @importFrom graphics lines plot.default text
#' @importFrom utils combn write.table read.csv write.csv
#' @importFrom stats runif

#' @name TreeMan-class
#' @title TreeMan-class
#' @aliases TreeMan-method
#' @description S4 class for representing phylogenetic trees as a list of nodes.
#' @param x \code{TreeMan} object
#' @param i node ID or slot name
#' @param object \code{TreeMan} object
#' @param max.level \code{str()} maximum number of levels to show
#' @param ... additional tree objects
#' @param j missing
#' @param drop missing
#' @slot ndlst list of nodes
#' @slot nds vector of node ids that are internal nodes
#' @slot nnds numeric of number of internal nodes in tree
#' @slot tips vector of node ids that are tips
#' @slot ntips numeric of number of internal nodes in tree
#' @slot all vector of all node ids
#' @slot nall numeric of number of all nodes in tree
#' @slot pd numeric of total branch length of tree
#' @slot tinds indexes of all tip nodes in tree
#' @slot prinds indexes of all pre-nodes in tree
#' @slot wspn logical, do nodes have spans
#' @slot wtxnyms logical, do nodes have txnyms
#' @slot ply logical, is tree bifurcating
#' @slot root character of node id of root, if no root then empty character
#' @slot updtd logical, if tree slots have been updated since initiation or change
#' @slot othr_slt_nms vector, character list of additional data slots added to nodes 
#' @slot ndmtrx matrix, T/Fs representing tree structure
#' @details
#' A \code{TreeMan} object holds a list of nodes. The idea of the \code{TreeMan}
#' class is to make adding and removing nodes as similar as possible to adding
#' and removing elements in a list. Note that internal nodes and tips are
#' both considered nodes. Trees can be polytomous but not unrooted.
#' 
#' 
#' Each node within the \code{TreeMan} \code{ndlst} contains the following data slots:
#' \itemize{
#'    \item \code{id}: character string for the node ID
#'    \item \code{txnym}: name of taxonomic clade (optional)
#'    \item \code{spn}: length of the preceding branch
#'    \item \code{prid}: ID of the immediately preceding node, NULL if root
#'    \item \code{ptid}: IDs of the immediately connecting nodes
#' }
#' 
#' See below in 'Examples' for these methods in use.
#' @seealso
#' \code{\link{randTree}}, \code{\link{Node-class}},
#' \code{\link{phylo-to-TreeMan}}, \code{\link{TreeMan-to-phylo}}
#' @examples
#' library(treeman)
#' # Generate random tree
#' tree <- randTree(10)
#' # Print to get basic stats
#' summary(tree)
#' # Slots....
#' tree['tips']   # return all tips IDs
#' tree['nds']    # return all internal node IDs
#' tree['ntips']  # count all tips
#' tree['nnds']   # count all internal nodes
#' tree['root']   # identify root node
#' tree[['t1']]   # return t1 node object
#' tree['pd']     # return phylogenetic diversity
#' tree['ply']    # is polytomous?
#' # Additional special slots
#' tree['age']   # get tree's age
#' tree['ultr']  # determine if tree is ultrametric
#' # Because all nodes are lists with metadata we can readily
#' #  get specific information on nodes of interest
#' nd <- tree[['n2']]
#' summary(nd)
#' # And then use the same syntax for the tree
#' nd['nkids']  # .... nkids, pd, etc.
#' 
#' # Convert to phylo and plot
#' library(ape)
#' tree <- as(tree, 'phylo')
#' plot(tree)
#' @exportClass TreeMan
setClass('TreeMan', representation=representation(
  ndlst='list',          # list of node lists
  nds='vector',          # vector of node ids that are internal nodes
  nnds='numeric',        # numeric of number of internal nodes in tree
  tips='vector',         # vector of node ids that are tips
  ntips='numeric',       # numeric of number of internal nodes in tree
  all='vector',          # vector of all Node ids
  nall='numeric',        # numeric of number of all nodes in tree
  pd='numeric',          # numeric of total branch length of tree
  wspn='logical',        # logical, do all nodes have spans
  wtxnyms='logical',     # logical, do nodes txnyms
  ply='logical',         # logical, is tree bifurcating
  updtd='logical',       # logical, if tree slots has been updated since a change
  ndmtrx='ANY',          # bigmemory matrix of logicals
  tinds='vector',        # indexes of tip nodes
  prinds='vector',       # indexes of pre-nodes
  root='character',      # character of node id of root, if no root then empty character
  othr_slt_nms='vector'),# if new slots added to node, list them here
  validity=fastCheckTreeMan)

# Accessor methods
#' @rdname TreeMan-class
#' @exportMethod [[
setMethod('[[', c('TreeMan', 'character'),
          function(x, i) {
            if(!i %in% names(x@ndlst)) {
              srch_trm <- gsub(' ', '_', i)  # usual mistake
              pssbls <- which(agrepl(srch_trm, names(x@ndlst), ignore.case=TRUE,
                                     max.distance=0.25))
              pssbls <- names(x@ndlst)[pssbls]
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
            .newNd(x, i)
          })
#' @rdname TreeMan-class
setMethod('[', c('TreeMan', 'character', 'missing', 'missing'),
          function(x, i, j, ..., drop=TRUE) {
            slt_nms <- slotNames(x)
            slt_nms <- slt_nms[slt_nms != 'ndlst']
            slt_nms <- slt_nms[slt_nms != 'ndmtrx']
            slt_nms <- slt_nms[slt_nms != 'tinds']
            slt_nms <- slt_nms[slt_nms != 'prinds']
            # ultr is special, shouldn't be updated when updateSlts()
            # too slow to calculate. Instead only calc if called.
            if(i == 'ultr') {
              return(isUltrmtrc(x))
            }
            if(i == 'age') {
              return(getAge(x))
            }
            if(!i %in% slt_nms) {
              slt_nms <- paste0(c(slt_nms, 'ultr', 'age'), collapse=', ')
              stop(paste0('`', i, '` not a tree slot. Available slots: ',
                          slt_nms))
            }
            slot(x, i)
          })

# display methods
#' @rdname TreeMan-class
#' @exportMethod as.character
setMethod('as.character', c('x'='TreeMan'),
          function(x) {
            paste0('TreeMan Object of [', length(x@tips),'] tips')
          })
#' @rdname TreeMan-class
#' @exportMethod show
setMethod('show', 'TreeMan',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname TreeMan-class
#' @exportMethod print
setMethod('print', 'TreeMan',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname TreeMan-class
#' @exportMethod str
setMethod('str', c('object'='TreeMan'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
#' @rdname TreeMan-class
#' @exportMethod summary
setMethod('summary', c('object'='TreeMan'),
          function(object){
            if(!fastCheckTreeMan(object)) {
              stop("Tree is corrupted. Run `checkNdlst()` to see how.")
            }
            if(!object@updtd) {
              stop("Tree is not updated since change or initiation. Use `updateSlts()`")
            }
            msg <- 'Tree (TreeMan Object):\n'
            msg <- paste0(msg, '  + ', object@ntips, ' tips\n')
            msg <- paste0(msg, '  + ', object@nnds, ' internal nodes\n')
            if(!is.null(object@ndmtrx)) {
              msg <- paste0(msg, '  + With node matrix\n')
            }
            if(object@wtxnyms) {
              msg <- paste0(msg, '  + With taxonomic names\n')
            }
            if(object@ply) {
              msg <- paste0(msg, '  + Polytomous\n')
            } else {
              msg <- paste0(msg, '  + Binary\n')
            }
            if(length(object@root) == 0) {
              if(!object@wspn) {
                msg <- paste0(msg, '  + Unrooted and without node spans\n')
              } else {
                msg <- paste0(msg, '  + Unrooted, with node spans\n')
                msg <- paste0(msg, '  + PD ', signif(object@pd, 3), '\n')
              }
            } else {
              if(object@wspn) {
                msg <- paste0(msg, '  + PD ', signif(object@pd, 3), '\n')
              } else {
                msg <- paste0(msg, '  + Without node spans\n')
              }
              msg <- paste0(msg, '  + Root node is \"', object@root, '\"\n')
            }
            if(length(object@othr_slt_nms) > 0) {
              msg <- paste0(msg, '  + With additional node slots:\n')
              for(slt_nm in object@othr_slt_nms) {
                msg <- paste0(msg, '    [', slt_nm, ']\n')
              }
            }
            cat(msg)
          })