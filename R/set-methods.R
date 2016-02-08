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
  diff <- tree@nodelist[[id]][['span']] - val
  tree@nodelist[[id]][['span']] <- diff
  tree@nodelist <- .updateNode(tree@nodelist, id, tree@root)
  tree@nodelist[[id]][['span']] <- val
  .updateSlots(tree)
}

setNodesSpan <- function(tree, ids, vals, ...) {
  .setNodeSpan <- function(id, val) {
    diff <- ndlst[[id]][['span']] - val
    ndlst[[id]][['span']] <- diff
    ndlst <- .updateNode(ndlst, id, tree@root)
    ndlst[[id]][['span']] <- val
    ndlst <<- ndlst
  }
  .nullify <- function(nd) {
    nd[['span']] <- NULL
    nd[['pd']] <- NULL
    nd[['prdst']] <- NULL
    nd
  }
  ndlst <- tree@nodelist
  if(is.null(vals)) {
    ndlst <- llply(ndlst[ids], .fun=.nullify)
  } else {
    m_ply(ids, .fun=.setNodeSpan, val=vals, ...)
  }
  tree@nodelist <- ndlst
  .updateSlots(tree)
}

setTol <- function(tree, tol) {
  tree@tol <- tol
  .updateSlots(tree)
}

setNodeID <- function(tree, id, val) {
  setNodesID(tree, id, val)
}

setNodesID <- function(tree, ids, vals, ...) {
  .rplcS4 <- function(slt) {
    if(any(slot(tree, slt) %in% ids)) {
      mtchs <- match(slot(tree, slt), ids)
      return(vals[mtchs])
    } else {
      return(slot(tree, slt))
    }
  }
  .run <-function(i) {
    .rplc <- function(slt) {
      if(any(nd[[slt]] %in% ids)) {
        mtchs <- match(nd[[slt]], ids)
        nd[[slt]] <- vals[mtchs]
      }
      nd
    }
    nd <- tree@nodelist[[i]]
    nd <- .rplc("id")
    nd <- .rplc("ptid")
    nd <- .rplc("prid")
    nd <- .rplc("kids")
    tree@nodelist[[i]] <<- nd
    NULL
  }
  l_data <- data.frame(i=1:length(tree@nodelist))
  m_ply(.data=l_data, .fun=.run, ...)
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
  .set <- function(i, val) {
    nd[[i]][['other']] <- val
    nd
  }
  l_data <- data.frame(i=1:length(vals), val=vals,
                       stringAsFactors=FALSE)
  tree@nodelist <- mlply(.data=l_data, .fun=.set)
  tree
}
