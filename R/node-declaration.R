.newNode <- function(tree, node) {
  node <- tree@nodelist[[node]]
  if(is.null(node[['span']]) | !tree@wspn) {
    span <- numeric()
  } else {
    span <- node[['span']]
  }
  if(length(tree@age) > 0) {
    age <- tree@age - node[['predist']]
  } else {
    age <- numeric()
  }
  new('Node', id=node[['id']], span=span, pre=as.character(node[['pre']]),
     post=as.character(node[['post']]), children=as.character(node[['children']]),
     nchildren=length(as.character(node[['children']])), pd=node[['pd']],
     predist=node[['predist']], root=tree@root == node[['id']],
     age=age, tip=length(node[['post']]) == 0)
}

setClass ('Node', representation=representation (
  id='character',        # unique ID for node in tree@nodelist
  span='numeric',        # length of preceding branch
  pre='character',       # parent node ID
  post='vector',         # child node IDs
  children='vector',     # descending tip IDs
  nchildren='numeric',   # number of descending tips
  pd='numeric',          # total branch length represented by node
  predist='numeric',     # total branch length of connected pres
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
             if(!x@root) {
               msg <- paste0(msg, '  + pre: \"', x@pre, '\"\n')
             }
             if(!x@tip) {
               msg <- paste0(msg, '  + post: \"', paste0(x@post, collapse='\", \"'), '\"\n')
               msg <- paste0(msg, '  + nchildren: ', length(x@children), '\n')
             }
             if(length(x@span) > 0) {
               if(!x@root) {
                 msg <- paste0(msg, '  + span: ', signif(x@span, 2), '\n')
               }
               if(length(x@age) > 0) {
                 msg <- paste0(msg, '  + age: ', signif(x@age, 2), '\n')
               } else {
                 msg <- paste0(msg, '  + predist: ', signif(x@predist, 2), '\n') 
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