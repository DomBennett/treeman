# #TODO: set_Taxonym, setRoot

setPD <- function(tree, val) {
  spans <- getNodesSlot(tree, ids=tree@all, name="span")
  spans <- spans/(tree@pd/val)
  tree <- setNodesSpan(tree, ids=tree@all, vals=spans)
  tree
}

setAge <- function(tree, val) {
  spans <- getNodesSlot(tree, ids=tree@all, name="span")
  spans <- spans/(tree@age/val)
  tree <- setNodesSpan(tree, ids=tree@all, vals=spans)
  tree
}

setNodeSpan <- function(tree, id, val) {
  .updateSlots(.setNodeSpan(tree, id, val))
}

# If vals=NULL, all node spans are set to NULL
# TODO: separate setNodeSpan and setNodesSpan, more efficient
setNodesSpan <- function(tree, ids, vals) {
  .run <- function(i) {
    tree <<- .setNodeSpan(tree, id=ids[i], val=vals[i])
    NULL
  }
  .nullify <- function(id) {
    tree@nodelist[[id]][['span']] <<- NULL
    tree@nodelist[[id]][['pd']] <<- NULL
    tree@nodelist[[id]][['prdst']] <<- NULL
  }
  if(is.null(vals)) {
    sapply(names(tree@nodelist), .nullify)
  } else {
    sapply(1:length(ids), .run)
  }
  .updateSlots(tree)
}

.setNodeSpan <- function(tree, id, val) {
  .pd <- function(prid, tree) {
    tree@nodelist[[prid]][['pd']] <- tree@nodelist[[prid]][['pd']] -
      pspn + val
    if(!is.null(tree@nodelist[[prid]][['prid']])) {
      tree <- .pd(tree@nodelist[[prid]][['prid']], tree)
    }
    tree
  }
  .prdst <- function(ptid, tree) {
    tree@nodelist[[ptid]][['prdst']] <- tree@nodelist[[ptid]][['prdst']] -
      pspn + val
    if(!is.null(tree@nodelist[[ptid]][['ptid']])) {
      for(ptid in tree@nodelist[[ptid]][['ptid']]) {
        tree <- .prdst(ptid, tree)
      }
    }
    tree
  }
  pspn <- tree@nodelist[[id]][['span']]
  tree@nodelist[[id]][['span']] <- val
  if(!is.null(tree@nodelist[[id]][['prid']])) {
    tree <- .pd(tree@nodelist[[id]][['prid']], tree)
  }
  tree <- .prdst(id, tree)
  tree
}

setTol <- function(tree, tol) {
  tree@tol <- tol
  .updateSlots(tree)
}

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


setNodeOther <- function(tree, id, val) {
  tree@nodelist[[id]]['other'] <- val
  tree
}

setNodesOther <- function(tree, ids, vals) {
  .set <- function(i) {
    tree <<- setNode(tree, ids[i], vals[i])
    NULL
  }
  sapply(1:length(ids), .set)
}
