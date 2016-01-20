# #TODO: set_ID, set_Span, set_Taxonym, set_Span, setTol
# .rename <- function(x, old, nw) {
#   .run <-function(nd) {
#     .rplc <- function(slt) {
#       if(any(nd[[slt]] %in% old)) {
#         mtchs <- match(nd[[slt]], old)
#         nd[[slt]] <- nw[mtchs]
#       }
#       nd
#     }
#     nd <- .rplc("id")
#     nd <- .rplc("postnode")
#     nd <- .rplc("prenode")
#     nd <- .rplc("children")
#     nd
#   }
#   x@nodelist <- lapply(x@nodelist, .run)
#   x
# }
# #TODO: move this to set-methods
# setGeneric("tips<-", signature=c("x"),
#            function(x, value) {
#              standardGeneric("tips<-")
#            })
# setReplaceMethod("tips", "TreeMan",
#                  function(x, value) {
#                    if(any(duplicated(value))) {
#                      stop('Tip names must be unique')
#                    }
#                    old_tips <- x@tips
#                    n <- length(old_tips)
#                    if(n != length(value)) {
#                      stop('Incorrect number of replacement tips')
#                    }
#                    x <- .rename(x, old_tips, value)
#                    mis <- match(old_tips, names(x@nodelist))
#                    names(x@nodelist)[mis] <- value
#                    x <- .update(x)
#                  })
# setGeneric("nodes<-", signature=c("x"),
#            function(x, value) {
#              standardGeneric("nodes<-")
#            })
# setReplaceMethod("nodes", "TreeMan",
#                  function(x, value) {
#                    if(any(duplicated(value))) {
#                      stop('Node names must be unique')
#                    }
#                    old_nodes <- x@nodes
#                    x <- .rename(x, old_nodes, value)
#                    x@root <- value[x@root == old_nodes]
#                    mis <- match(old_nodes, names(x@nodelist))
#                    names(x@nodelist)[mis] <- value
#                    .update(x)
#                  })
# 
# Accessor method
setGeneric('setTol', signature=c('x', 'n'),
           function(x, n) {
             genericFunction('setTol')
           })
setMethod('setTol', c('TreeMan', 'numeric'),
          function(x, n){
            x@tol <- n
            .update(x)
          })

setNode <- function(tree, id, name, value) {
  tree@nodelist[[id]][name] <- value
  tree
}

setNodes <- function(tree, ids, name, values) {
  .set <- function(i) {
    tree <<- setNode(tree, ids[i], name, values[i])
    NULL
  }
  sapply(1:length(ids), .set)
  treeman:::.update(tree)
}
