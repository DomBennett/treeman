# #TODO: set_Span, set_Taxonym, setTol

setGeneric('setTol', signature=c('tree', 'tol'),
           function(tree, tol) {
             genericFunction('setTol')
           })
setMethod('setTol', c('TreeMan', 'numeric'),
          function(tree, tol){
            tree@tol <- tol
            .update(tree)
          })

# setRoot
# setAge
# setPD
# setNodeTaxonym

setNodeID <- function(tree, id, val) {
  setNodesID(tree, id, val)
}

setNodesID <- function(tree, ids, vals) {
  .rplcS4 <- function(slt) {
    if(any(slot(tree, slt) %in% ids)) {
      mtchs <- match(slot(tree, slt), ids)
      return(vals[mtchs])
    } else {
      return(slot(tree, slt))
    }
  }
  .run <-function(nd) {
    .rplc <- function(slt) {
      if(any(nd[[slt]] %in% ids)) {
        mtchs <- match(nd[[slt]], ids)
        nd[[slt]] <- vals[mtchs]
      }
      nd
    }
    nd <- .rplc("id")
    nd <- .rplc("ptid")
    nd <- .rplc("prid")
    nd <- .rplc("children")
    nd
  }
  tree@nodelist <- lapply(tree@nodelist, .run)
  tree@tips <- .rplcS4('tips')
  tree@nodes <- .rplcS4('nodes')
  tree@ext <- .rplcS4('ext')
  tree@exc <- .rplcS4('exc')
  tree@root <- .rplcS4('root')
  tree
}


setNodeOther <- function(tree, id, value) {
  tree@nodelist[[id]]['other'] <- value
  tree
}

setNodesOther <- function(tree, ids, values) {
  .set <- function(i) {
    tree <<- setNode(tree, ids[i], values[i])
    NULL
  }
  sapply(1:length(ids), .set)
}
