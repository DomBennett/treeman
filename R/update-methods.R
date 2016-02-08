# Convert all from mlply to llply, much faster
.dwndateNode <- function(ndlst, nid, rid) {
  .add <- function(id) {
    ndlst[[id]][['pd']] <<- ndlst[[id]][['pd']] - nd_span
  }
  nd_span <- ndlst[[nid]][['span']]
  prids <- ndlst[[nid]][['prid']]
  mlply(prids, .fun=.add)
  ndlst
}

.updateNode <- function(ndlst, nid, rid) {
  .add <- function(id) {
    ndlst[[id]][['pd']] <<- ndlst[[id]][['pd']] + nd_span
    prdst <<- ndlst[[id]][['span']] + prdst
  }
  nd_span <- ndlst[[nid]][['span']]
  prdst <- nd_span
  prids <- ndlst[[nid]][['prid']]
  mlply(prids, .fun=.add)
  ndlst[[nid]][['prdst']] <- prdst
  ndlst
}

.dwndateTip <- function(ndlst, tid, rid) {
  .add <- function(id) {
    kids <- ndlst[[id]][['kids']]
    ndlst[[id]][['kids']] <<- kids[kids != tid]
    ndlst[[id]][['pd']] <<- ndlst[[id]][['pd']] - tp_span
  }
  tp_span <- ndlst[[tid]][['span']]
  prdst <- 0
  prids <- ndlst[[tid]][['prid']]
  mlply(prids, .fun=.add)
  ndlst[[tid]][['prdst']] <- prdst
  ndlst
}

.updateTip <- function(ndlst, tid, rid) {
  .add <- function(nd) {
    kids <- nd[['kids']]
    nd[['kids']] <- c(kids, tid)
    nd[['pd']] <- nd[['pd']] + tp_span
    prdst <<- nd[['span']] + prdst
    nd
  }
  tp_span <- ndlst[[tid]][['span']]
  prids <- ndlst[[tid]][['prid']]
  prdst <- tp_span
  ndlst[prids] <- llply(ndlst[prids], .add)
  ndlst[[tid]][['prdst']] <- prdst
  ndlst
}

.globalUpdateAll <- function(ndlst) {
  tip <- function(tid) {
    ndlst <- .updateTip(ndlst, tid, rid)
    ndlst <<- ndlst
  }
  node <- function(nid) {
    ndlst <- .updateNode(ndlst, nid, rid)
    ndlst <<- ndlst
  }
  wo_pstnds <- sapply(ndlst, function(n) length(n[['ptid']]) == 0)
  wo_prnds <- sapply(ndlst, function(n) length(n[['prid']]) == 0)
  nids <- names(ndlst)[(!wo_pstnds) & (!wo_prnds)]
  tids <- names(ndlst)[wo_pstnds]
  rid <- names(ndlst)[wo_prnds]
  l_data <- data.frame(tid=tids, stringsAsFactors=FALSE)
  m_ply(.data=l_data, .fun=tip)
  l_data <- data.frame(nid=nids, stringsAsFactors=FALSE)
  m_ply(.data=l_data, .fun=node)
  ndlst
}

.updateKids <- function(ndlst, tid, rid) {
  .add <- function(id) {
    kids <- ndlst[[id]][['kids']]
    ndlst[[id]][['kids']] <<- c(kids, tid)
  }
  prids <- ndlst[[tid]][['prid']]
  mlply(prids, .fun=.add)
  ndlst
}

.dwndateKids <- function(ndlst, tid, rid) {
  .add <- function(id) {
    kids <- ndlst[[id]][['kids']]
    ndlst[[id]][['kids']] <<- kids[kids != tid]
  }
  prids <- ndlst[[tid]][['prid']]
  mlply(prids, .fun=.add)
  ndlst
}

.globalUpdateKids <- function(ndlst) {
  tip <- function(tid) {
    ndlst <- .updateKids(ndlst, tid, rid)
    ndlst <<- ndlst
  }
  wo_pstnds <- sapply(ndlst, function(n) length(n[['ptid']]) == 0)
  w_prnds <- sapply(ndlst, function(n) length(n[['prid']]) == 0)
  tids <- names(ndlst)[wo_pstnds]
  rid <- names(ndlst)[w_prnds]
  l_data <- data.frame(tid=tids, stringsAsFactors=FALSE)
  m_ply(.data=l_data, .fun=tip)
  ndlst
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