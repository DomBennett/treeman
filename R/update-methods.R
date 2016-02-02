.updateNodes <- function(tree_env, nid, rid) {
  .add <- function(tree_env) {
    id <- ndlst[[tree_env[['id']]]][['prid']]
    span <- tree_env[['ndlst']][[id]][['span']]
    tree_env[['prdst']] <- tree_env[['prdst']] + span
    pd <- tree_env[['ndlst']][[id]][['pd']] + nd_span
    tree_env[['ndlst']][[id]][['pd']] <- pd
    if(tree_env[['id']] != rid) {
      tree_env[['id']] <- id
      .add(tree_env)
    }
  }
  tree_env[['id']] <- tree_env[['ndlst']][[nid]][['prid']]
  tip_span <- tree_env[['ndlst']][[nid]][['span']]
  tree_env[['prdst']] <- nd_span
  finish <- try(stop(), silent=TRUE)
  while(is(finish, 'try-error')) {
    finish <- try(expr={
      .add(tree_env)
    }, silent=TRUE)
  }
  tree_env[['ndlst']][[nid]][['prdst']] <- tree_env[['prdst']]
  NULL
}

.updateTips <- function(tree_env, tid, rid) {
  .add <- function(tree_env) {
    id <- ndlst[[tree_env[['id']]]][['prid']]
    kids <- c(tree_env[['ndlst']][[id]][['children']], tid)
    tree_env[['ndlst']][[id]][['children']] <- kids
    span <- tree_env[['ndlst']][[id]][['span']]
    tree_env[['prdst']] <- tree_env[['prdst']] + span
    pd <- tree_env[['ndlst']][[id]][['pd']] + tp_span
    tree_env[['ndlst']][[id]][['pd']] <- pd
    if(tree_env[['id']] != rid) {
      tree_env[['id']] <- id
      .add(tree_env)
    }
  }
  tree_env[['id']] <- tree_env[['ndlst']][[tid]][['prid']]
  tp_span <- tree_env[['ndlst']][[tid]][['span']]
  tree_env[['prdst']] <- tp_span
  finish <- try(stop(), silent=TRUE)
  while(is(finish, 'try-error')) {
    finish <- try(expr={
      .add(tree_env)
    }, silent=TRUE)
  }
  tree_env[['ndlst']][[tid]][['prdst']] <- tree_env[['prdst']]
  NULL
}

.globalUpdateAll <- function(ndlst, tids, nids, rid, ...) {
  tree_env <- new.env()
  tree_env$ndlst <- ndlst
  l_data <- data.frame(tid=tids, stringsAsFactors=FALSE)
  m_ply(.data=l_data, .fun=.updateTips, tree_env=tree_env, rid=rid, ...)
  l_data <- data.frame(nid=nids, stringsAsFactors=FALSE)
  m_ply(.data=l_data, .fun=.updateNodes, tree_env=tree_env, rid=rid, ...)
  tree_env$ndlst
}

.localUpdateTip <- function(ndlst, tid, rid) {
  # update all node slots for a tip
  tree_env <- new.env()
  tree_env$ndlst <- ndlst
  .updateTip(tree_env, tid, rid)
  tree_env$ndlst
}

.localUpdateNode <- function(ndlst, nid, rid) {
  # update all node slots for an internal node
  tree_env <- new.env()
  tree_env$ndlst <- ndlst
  .updateNode(tree_env, nid, rid)
  tree_env$ndlst
}

.updateChildren <- function(tree_env, tid, rid) {
  .add <- function(tree_env) {
    id <- tree_env[['ndlst']][[tree_env[['id']]]][['prid']]
    kids <- c(tree_env[['ndlst']][[id]][['children']], tid)
    tree_env[['ndlst']][[id]][['children']] <- kids
    if(tree_env[['id']] != rid) {
      tree_env[['id']] <- id
      .add(tree_env)
    }
  }
  tree_env[['id']] <- tree_env[['ndlst']][[tid]][['prid']]
  finish <- try(stop(), silent=TRUE)
  while(is(finish, 'try-error')) {
    finish <- try(expr={
      .add(tree_env)
    }, silent=TRUE)
  }
  NULL
}

.globalUpdateChildren <- function(ndlst, tids, rid, ...) {
  # add children to all nodes in tree
  tree_env <- new.env()
  tree_env$ndlst <- ndlst
  l_data <- data.frame(tid=tids, stringsAsFactors=FALSE)
  m_ply(.data=l_data, .fun=.updateChildren, tree_env=tree_env, rid=rid, ...)
  tree_env$ndlst
}

.localUpdateChildren <- function(ndlst, tid, rid) {
  # run this to add just children slot for new tip
  tree_env <- new.env()
  tree_env$ndlst <- ndlst
  .updateChildren(tree_env, tid, rid)
  tree_env$ndlst
}

.updateSlots <- function(tree) {
  wo_pstndes <- sapply(tree@nodelist,
                       function(n) length(n[['ptid']]) == 0)
  tree@tips <- sort(names(wo_pstndes)[wo_pstndes])
  tree@ntips <- length(tree@tips)
  tree@nodes <- sort(names(wo_pstndes)[!wo_pstndes])
  tree@nnodes <- length(tree@nodes)
  tree@all <- c(tree@tips, tree@nodes)
  tree@nall <- length(tree@all)
  wspn <- names(tree@nodelist)[names(tree@nodelist) != tree@root]
  tree@wspn <- all(sapply(tree@nodelist[wspn], function(n) !is.null(n[['span']])))
  if(tree@wspn) {
    if(!is.null(tree@root)) {
      tree@age <- max(sapply(tree@nodelist[wspn], function(n) n[['prdst']]))
      extant_is <- unlist(sapply(tree@tips, function(i) {
        (tree@age - tree@nodelist[[i]][['prdst']]) <= tree@tol}))
      tree@ext <- names(extant_is)[extant_is]
      tree@exc <- tree@tips[!tree@tips %in% tree@ext]
      tree@ultr <- all(tree@tips %in% tree@ext)
    }
    tree@pd <- tree@nodelist[[tree@root]][['pd']]
  } else {
    tree@age <- tree@pd <- numeric()
    tree@ext <- tree@ext <- vector()
    tree@ultr <- logical()
  }
  tree@ply <- any(sapply(tree@nodelist, function(n) length(n[['ptid']]) > 2))
  initialize(tree)
}