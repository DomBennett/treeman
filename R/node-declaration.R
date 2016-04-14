.newNode <- function(tree, node) {
  node <- tree@nodelist[[node]]
  if(is.null(node[['span']]) | !tree@wspn) {
    span <- pd <- prdst <- numeric()
  } else {
    span <- node[['span']]
    pd <- node[['pd']]
    prdst <- node[['prdst']]
  }
  if(length(tree@age) > 0) {
    age <- tree@age - node[['prdst']]
  } else {
    age <- numeric()
  }
  if(is.null(node[['txnym']])) {
    txnym <- vector()
  } else {
    txnym <- node[['txnym']]
  }
  new('Node', id=node[['id']], span=span, prid=as.character(node[['prid']][1]),
     ptid=as.character(node[['ptid']]), kids=as.character(node[['kids']]),
     nkids=length(as.character(node[['kids']])), pd=pd, txnym=txnym,
     prdst=prdst, root=tree@root == node[['id']],
     age=age, tip=length(node[['ptid']]) == 0)
}

#' @name Node
#' @title S4 class for displaying nodes
#' @description The \code{Node} class is used to display node information.
#' It is only generated when a user implements the \code{[[]]} on a tree.
#' @slot id unique ID for node in tree['nodelist']
#' @slot span length of preceding branch
#' @slot prid parent node ID
#' @slot ptid child node ID
#' @slot kids descending tip IDs
#' @slot nkids number of descending tip IDs
#' @slot txnym list of associated taxonyms
#' @slot pd total branch length represented by node
#' @slot prdst total branch length of connected prids
#' @slot age age of node in tree
#' @slot root T/F root node?
#' @slot tip T/F tip node?
#' @exportClass Node
#' @seealso 
#' \code{\link{cTrees}}
setClass ('Node', representation=representation (
  id='character',        # unique ID for node in tree@nodelist
  span='numeric',        # length of preceding branch
  prid='character',      # parent node ID
  ptid='vector',         # child node IDs
  kids='vector',         # descending tip IDs
  nkids='numeric',       # number of descending tips
  txnym="vector",        # list of associated taxonyms
  pd='numeric',          # total branch length represented by node
  prdst='numeric',       # total branch length of connected pres
  age='numeric',         # age of node in tree
  root='logical',        # T/F root node?
  tip='logical')         # T/F tip node?
)

setMethod ('as.character', c('x'='Node'),
           function(x) {
             x@id
           })
setMethod ('show', 'Node',
           function(object){
             print (object)
           })
setGeneric ('print')
setMethod ('print', c('x'='Node'),
           function(x){
             if(x@root) {
               msg <- paste0('Node (root node):\n')
             } else if (x@tip){
               msg <- paste0('Node (tip node):\n')
             } else {
               msg <- paste0('Node (internal node):\n')
             }
             msg <- paste0(msg, '  + ID: \"', x@id, '\"\n')
             if(length(x@txnym) > 0) {
               msg <- paste0(msg, '  + txnym: \"', paste0(x@txnym, collapse='\", \"'), '\"\n')
             }
             if(!x@root) {
               msg <- paste0(msg, '  + prid: \"', x@prid, '\"\n')
             }
             if(!x@tip) {
               msg <- paste0(msg, '  + ptid: \"', paste0(x@ptid, collapse='\", \"'), '\"\n')
               msg <- paste0(msg, '  + nkids: ', length(x@kids), '\n')
             }
             if(length(x@span) > 0) {
               if(!x@root) {
                 msg <- paste0(msg, '  + span: ', signif(x@span, 2), '\n')
               }
               if(length(x@age) > 0) {
                 msg <- paste0(msg, '  + age: ', signif(x@age, 2), '\n')
               } else {
                 msg <- paste0(msg, '  + predist: ', signif(x@prdst, 2), '\n') 
               }
               msg <- paste0(msg, '  + pd: ', signif(x@pd, 2), '\n')
             }
             cat (msg)
           })
setMethod('[', c('Node', 'character'),
          function(x, i) {
            if(!i %in% slotNames(x)) {
              stop(paste0(i, '  not in node'))
            }
            slot(x, i)
          })