.newNode <- function(tree, node) {
  node <- tree@nodelist[[node]]
  new('Node', id=node$id, span=node$span, prenode=as.character(node$prenode),
     postnode=as.character(node$postnode), children=as.character(node$children),
     pd=node$pd, predist=node$predist, tree_age=tree@age, root=tree@root == node$id,
     tip=length(node$postnode) == 0)
}

setClass ('Node', representation=representation (
  id='character',        # unique ID for node in tree@nodelist
  span='numeric',        # length of preceding branch
  prenode='character',   # parent node ID
  postnode='vector',     # child node IDs
  children='vector',     # descending tip IDs
  pd='numeric',          # total branch length represented by node
  predist='numeric',     # total branch length of connected prenodes
  tree_age='ANY',        # age of original tree
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
             msg <- paste0(msg, '-- ID: [`', x@id, '`]\n')
             if(!x@root) {
               msg <- paste0(msg, '-- prenode: [`', x@prenode, '`]\n')
             }
             if(!x@tip) {
               msg <- paste0(msg, '-- postnode: [`', paste0(x@postnode, collapse='`, `'), '`]\n')
               msg <- paste0(msg, '-- nChildren: [', length(x@children), ']\n')
             }
             if(!x@root) {
               msg <- paste0(msg, '-- span: [', signif(x@span, 2), ']\n')
             }
             if(is.numeric(x@tree_age)) {
               msg <- paste0(msg, '-- age: [', signif(x@tree_age - x@predist, 2), ']\n')
             } else {
               msg <- paste0(msg, '-- predist: [', signif(x@predist, 2), ']\n') 
             }
             msg <- paste0(msg, '-- PD: [', signif(x@pd, 2), ']\n')
             cat (msg)
           })
