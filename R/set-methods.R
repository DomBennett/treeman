# #TODO
setNodeTaxonym <- function(...) {
  cat('Sorry not yet implemented!\n')
}

setNodesTaxonym <- function(...) {
  cat('Sorry not yet implemented!\n')
}

setRoot <- function(...) {
  cat('Sorry not yet implemented!\n')
}

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
  .ptnd <- function(nd) {
    nd[['prdst']] <- nd[['prdst']] + diff
    nd
  }
  # reset node using diff
  diff <- val - tree@nodelist[[id]][['span']]
  tree@nodelist[[id]][['span']] <- diff
  tree@nodelist <- .updateNode(tree@nodelist, id, tree@root)
  # adjust any pstnds
  ptids <- getNodePtid(tree, id=id)
  ptids <- ptids[-(length(ptids))]
  tree@nodelist[ptids] <- llply(tree@nodelist[ptids], .fun=.ptnd)
  # update nd
  tree@nodelist[[id]][['span']] <- val
  tree@nodelist[[id]][['prdst']] <- tree@nodelist[[id]][['prdst']] +
    val + abs(diff)
  .updateTreeSlots(tree)
}

setNodesSpan <- function(tree, ids, vals, ...) {
  .nullify <- function(nd) {
    nd[['span']] <- NULL
    nd[['pd']] <- NULL
    nd[['prdst']] <- NULL
    nd
  }
  .reset <- function(id, span) {
    ndlst[[id]][['span']] <- span
    ndlst[[id]][['pd']] <- 0
    ndlst[[id]][['prdst']] <- 0
    ndlst[[id]]
  }
  ndlst <- tree@nodelist
  if(is.null(vals)) {
    ndlst <- llply(ndlst[tree@all], .fun=.nullify, ...)
  } else {
    spans <- getNodesSlot(tree, name='span', ids=tree@all)
    spans[match(ids, tree@all)] <- vals
    l_data <- data.frame(id=tree@all, span=spans, stringsAsFactors=FALSE)
    ndlst <- mlply(l_data, .fun=.reset)
    ndlst <- ndlst[1:length(ndlst)]
    names(ndlst) <- tree@all
    ndlst <- .globalUpdateAll(ndlst, just_spn_data=TRUE)
  }
  tree@nodelist <- ndlst
  .updateTreeSlots(tree)
}

setTol <- function(tree, tol) {
  tree@tol <- tol
  .updateTreeSlots(tree)
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
