.newNd <- function(tree, id) {
  nd <- tree@ndlst[[id]]
  if(!tree@wspn) {
    spn <- pd <- prdst <- numeric()
  } else {
    spn <- nd[['spn']]
    pd <- .getNdPDFrmLst(tree@ndlst, prinds=tree@prinds,
                                   id=id)
    prdst <- .getNdPrdstsFrmLst(tree@ndlst, prinds=tree@prinds,
                                id=id)
  }
  if(is.null(nd[['txnym']])) {
    txnym <- vector()
  } else {
    txnym <- nd[['txnym']]
  }
  kids <- .getNdKidsFrmLst(tree@ndlst, prinds=tree@prinds, id=id, tinds=tree@tinds)
  new('Node', id=nd[['id']], spn=spn, prid=as.character(nd[['prid']][1]),
     ptid=as.character(nd[['ptid']]), kids=as.character(kids),
     nkids=length(kids), pd=pd, txnym=txnym, prdst=prdst,
     root=tree@root == nd[['id']], tip=length(nd[['ptid']]) == 0)
}

#' @name Node-class
#' @aliases Node-method
#' @param x \code{Node} object
#' @param object \code{Node} object
#' @param i slot name
#' @param j missing
#' @param ... missing
#' @param drop missing
#' @title Node-class
#' @description The \code{Node} is an S4 class used for displaying node information.
#' It is only generated when a user implements the \code{[[]]} on a tree. Information
#' is only accurate if tree has been updated with \code{updateTree()}.
#' @slot id unique ID for node in tree['ndlst']
#' @slot spn length of preceding branch
#' @slot prid parent node ID
#' @slot ptid child node ID
#' @slot kids descending tip IDs
#' @slot nkids number of descending tip IDs
#' @slot txnym list of associated taxonyms
#' @slot pd total branch length represented by node
#' @slot prdst total branch length of connected prids
#' @slot root T/F root node?
#' @slot tip T/F tip node?
#' @exportClass Node
#' @seealso 
#' \code{\link{TreeMan-class}}, \code{\link{TreeMen-class}}
setClass ('Node', representation=representation (
  id='character',        # unique ID for node in tree@nodelist
  spn='numeric',        # length of preceding branch
  prid='character',      # parent node ID
  ptid='vector',         # child node IDs
  kids='vector',         # descending tip IDs
  nkids='numeric',       # number of descending tips
  txnym="vector",        # list of associated taxonyms
  pd='numeric',          # total branch length represented by node
  prdst='numeric',       # total branch length of connected pres
  root='logical',        # T/F root node?
  tip='logical')         # T/F tip node?
)
#' @rdname Node-class
#' @exportMethod as.character
setMethod ('as.character', c('x'='Node'),
           function(x) {
             paste0('Node Obj. [ID=', x@id, ']')
           })
#' @rdname Node-class
#' @exportMethod show
setMethod ('show', 'Node',
           function(object){
             cat(summary(object))
           })
#' @rdname Node-class
#' @exportMethod print
setMethod ('print', 'Node',
           function(x){
             print(summary(x))
           })
#' @rdname Node-class
#' @exportMethod summary
setMethod ('summary', c('object'='Node'),
           function(object){
             if(object@root) {
               msg <- paste0('Node (root node):\n')
             } else if (object@tip){
               msg <- paste0('Node (tip node):\n')
             } else {
               msg <- paste0('Node (internal node):\n')
             }
             msg <- paste0(msg, '  + ID: \"', object@id, '\"\n')
             if(length(object@txnym) > 0) {
               msg <- paste0(msg, '  + txnym: \"', paste0(object@txnym, collapse='\", \"'), '\"\n')
             }
             if(!object@root) {
               msg <- paste0(msg, '  + prid: \"', object@prid, '\"\n')
             }
             if(!object@tip) {
               msg <- paste0(msg, '  + ptid: \"', paste0(object@ptid, collapse='\", \"'), '\"\n')
               msg <- paste0(msg, '  + nkids: ', length(object@kids), '\n')
             }
             if(length(object@spn) > 0) {
               if(!object@root) {
                 msg <- paste0(msg, '  + spn: ', signif(object@spn, 2), '\n')
               }
               msg <- paste0(msg, '  + predist: ', signif(object@prdst, 2), '\n') 
               msg <- paste0(msg, '  + pd: ', signif(object@pd, 2), '\n')
             }
             cat(msg)
           })

#' @rdname Node-class
#' @exportMethod [
setMethod('[', c('Node', 'character', 'missing', 'missing'),
          function(x, i, j, ..., drop=TRUE) {
            if(!i %in% slotNames(x)) {
              stop(paste0(i, '  not in node'))
            }
            slot(x, i)
          })