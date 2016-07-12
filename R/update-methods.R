
# TODO: need to rethink this, make it more logical
# -- update tips or nodes
# -- update upstream or downstream or all
# -- update prdst, pd, kids or all
# -- update with C

.getTreels <- function(ndlst) {
  # Breaks a ndlst into key tree elements
  # ids, tids (as integers), prids (as integers)
  # this is needed for running with C
  # get ids
  ids <- names(ndlst)
  # get prids
  prids <- sapply(ndlst, function(x) x[['prid']])
  prids <- match(prids, ids)
  prids[is.na(prids)] <- -1
  # get tids
  tids <- ids[sapply(ndlst, function(x) length(x[['ptid']]) == 0)]
  tids <- match(tids, ids)
  # get spns
  spns <- as.numeric(sapply(ndlst, function(x) x[['spn']]))
  spns[is.na(spns)] <- 0
  # return
  list('ids'=ids, 'tids'=tids, 'prids'=prids, 'spns'=spns)
}

#' @useDynLib treeman
#' @useDynLib treeman getKidsMat
#' @useDynLib treeman getPrdstVec
#' @useDynLib treeman getPdVec

.updateNdlst <- function(treels, spn_offset=0) {
  #### Complex function #####
  # (re)calculate tree slots from first principles
  # in case of any changes made or new tree to be created
  # takes tree elements (treels) and returns ndlst
  # assumes ids, spns, ptids and prids are correct
  # (re)calculates: kids, prdsts and pds, where available
  # takes spn_offset, the distance within tree
  ###########################
  # Internals
  .addwospn <- function(i) {
    nd <- vector("list", length=4)
    names(nd) <- c('id', 'ptid', 'prid', 'kids')
    nd[['id']] <- ids[i]
    nd[['kids']] <- tids[kids[i, ]]
    nd[['prid']] <- prids[i]
    ptids <- ids[prids == ids[i]]
    nd[['ptid']] <- ptids[!is.na(ptids)]
    nd
  }
  .addwspn <- function(i) {
    nd <- vector("list", length=7)
    names(nd) <- c('id', 'ptid', 'prid', 'kids',
                   'spn', 'pd', 'prdst')
    nd[['id']] <- ids[i]
    nd[['spn']] <- spns[i]
    nd[['prdst']] <- prdsts[i]
    nd[['pd']] <- pds[i]
    nd[['kids']] <- tids[kids[i, ]]
    nd[['prid']] <- prids[i]
    ptids <- ids[prids == ids[i]]
    nd[['ptid']] <- ptids[!is.na(ptids)]
    nd
  }
  # unpack
  ids <- treels[['ids']]
  nids <- as.integer(length(ids))
  tids <- as.integer(treels[['tids']])
  prids <- as.integer(treels[['prids']])
  spns <- as.numeric(treels[['spns']])
  # get kids
  kids <- .Call("getKidsMat", PACKAGE="treeman",
                nids, tids, prids)
  kids <- kids == 1
  if(sum(spns) > 0) {
    prdsts <- .Call("getPrdstVec", PACKAGE="treeman",
                    nids, prids, spns)
    prdsts <- prdsts + spn_offset
    pds <- .Call("getPdVec", PACKAGE="treeman",
                 nids, prids, spns)
  }
  # replace -1s with NAs
  prids[prids == -1] <- NA
  spns[spns == 0] <- NA
  # generate ndlst
  tids <- ids[tids]
  prids <- ids[prids]
  if(sum(is.na(spns)) > 1) {
    ndlst <- lapply(1:length(ids), .addwospn)
  } else {
    ndlst <- lapply(1:length(ids), .addwspn)
  }
  names(ndlst) <- ids
  ndlst
}

.updateNd <- function(tree, id) {
  # Update nodes from changed node id
  # get all pre-nodes
  prids <- getNdPrids(tree, id)
  # get all pstnds
  ptids <- getNdPtid(tree, id)
  # ids
  ids <- c(prids, id, ptids)
  # update
  tree@ndlst[ids] <-
    .updateNdlst(.getTreels(tree@ndlst[ids]))
  .updateTreeSlts(tree)
}

.updateNdsSlt <- function(ndlst, nids, updater) {
  # update nids using updater function
  ndlst[nids] <- plyr::llply(ndlst[nids], .fun=updater)
  ndlst
}

.dwndateNd <- function(ndlst, nid, rid) {
  .add <- function(nd) {
    nd[['pd']] <- nd[['pd']] - nd_spn
    nd[['kids']] <- nd[['kids']][nd[['kids']] != nid]
    nd
  }
  nd_spn <- ndlst[[nid]][['spn']]
  prids <- ndlst[[nid]][['prid']]
  ndlst[prids] <- plyr::llply(ndlst[prids], .fun=.add)
  ndlst
}

.updateNd <- function(ndlst, nid, rid) {
  .add <- function(nd) {
    nd[['pd']] <- nd[['pd']] + nd_spn
    prdst <<- nd[['spn']] + prdst
    nd
  }
  nd_spn <- ndlst[[nid]][['spn']]
  prdst <- nd_spn
  prids <- ndlst[[nid]][['prid']]
  ndlst[prids] <- plyr::llply(ndlst[prids], .add)
  ndlst[[nid]][['prdst']] <- prdst
  ndlst
}

.updateNd2 <- function(tree, nid) {
  .add <- function(nd) {
    nd[['pd']] <- nd[['pd']] + nd_spn
    prdst <<- nd[['spn']] + prdst
    nd
  }
  rid <- tree@root
  nd_spn <- tree@ndlst[[nid]][['spn']]
  prdst <- nd_spn
  prids <- getNdPrids(tree, nid)
  tree@ndlst[prids] <- plyr::llply(tree@ndlst[prids], .add)
  tree@ndlst[[nid]][['prdst']] <- prdst
  tree
}

.dwndateTip <- function(ndlst, tid, rid) {
  .add <- function(nd) {
    kids <- nd[['kids']]
    nd[['kids']] <- kids[kids != tid]
    nd[['pd']] <- nd[['pd']] - tp_spn
    nd
  }
  tp_spn <- ndlst[[tid]][['spn']]
  prids <- ndlst[[tid]][['prid']]
  ndlst[prids] <- plyr::llply(ndlst[prids], .fun=.add)
  ndlst
}

.updateTip <- function(ndlst, tid, rid) {
  .add <- function(nd) {
    kids <- nd[['kids']]
    nd[['kids']] <- c(kids, tid)
    nd[['pd']] <- nd[['pd']] + tp_spn
    prdst <<- nd[['spn']] + prdst
    nd
  }
  tp_spn <- ndlst[[tid]][['spn']]
  prids <- ndlst[[tid]][['prid']]
  prdst <- tp_spn
  ndlst[prids] <- plyr::llply(ndlst[prids], .add)
  ndlst[[tid]][['prdst']] <- prdst
  ndlst
}

.globalUpdateAll <- function(ndlst, just_spn_data=FALSE) {
  tip <- function(tid) {
    ndlst <- .updateTip(ndlst, tid, rid)
    ndlst <<- ndlst
  }
  nd <- function(nid) {
    ndlst <- .updateNd(ndlst, nid, rid)
    ndlst <<- ndlst
  }
  wo_prnds <- sapply(ndlst, function(n) length(n[['prid']]) == 0)
  if(!just_spn_data) {
    wo_pstnds <- sapply(ndlst, function(n) length(n[['ptid']]) == 0)
    nids <- names(ndlst)[(!wo_pstnds) & (!wo_prnds)]
    tids <- names(ndlst)[wo_pstnds]
    rid <- names(ndlst)[wo_prnds]
    l_data <- data.frame(tid=tids, stringsAsFactors=FALSE)
    plyr::m_ply(.data=l_data, .fun=tip)
  } else {
    # just run updateNd for all nodes if just spn data needs updating
    nids <- names(ndlst)[!wo_prnds]
    rid <- names(ndlst)[wo_prnds]
  }
  l_data <- data.frame(nid=nids, stringsAsFactors=FALSE)
  plyr::m_ply(.data=l_data, .fun=nd)
  ndlst
}

.updateKids <- function(ndlst, tid, rid) {
  .add <- function(nd) {
    kids <- nd[['kids']]
    nd[['kids']] <- c(kids, tid)
    nd
  }
  prids <- ndlst[[tid]][['prid']]
  ndlst[prids] <- plyr::llply(ndlst[prids], .fun=.add)
  ndlst
}

.dwndateKids <- function(ndlst, tid, rid) {
  .add <- function(nd) {
    kids <- nd[['kids']]
    nd[['kids']] <- kids[kids != tid]
    nd
  }
  prids <- ndlst[[tid]][['prid']]
  ndlst[prids] <- plyr::llply(ndlst[prids], .fun=.add)
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
  plyr::m_ply(.data=l_data, .fun=tip)
  ndlst
}

.updateTreeSlts <- function(tree) {
  wo_pstndes <- sapply(tree@ndlst,
                       function(n) length(n[['ptid']]) == 0)
  tree@tips <- sort(names(wo_pstndes)[wo_pstndes])
  tree@ntips <- length(tree@tips)
  tree@nds <- sort(names(wo_pstndes)[!wo_pstndes])
  tree@nnds <- length(tree@nds)
  tree@all <- c(tree@tips, tree@nds)
  tree@nall <- length(tree@all)
  if(length(tree@root) > 0) {
    wspn <- names(tree@ndlst)[names(tree@ndlst) != tree@root]
  } else {
    wspn <- names(tree@ndlst)
  }
  tree@wspn <- all(sapply(tree@ndlst[wspn], function(n) !is.null(n[['spn']])))
  if(tree@wspn) {
    if(length(tree@root) > 0) {
      tree@age <- max(sapply(tree@ndlst[wspn], function(n) n[['prdst']]))
      extant_is <- unlist(sapply(tree@tips, function(i) {
        (tree@age - tree@ndlst[[i]][['prdst']]) <= tree@tol}))
      tree@ext <- names(extant_is)[extant_is]
      tree@exc <- tree@tips[!tree@tips %in% tree@ext]
      tree@ultr <- all(tree@tips %in% tree@ext)
    } else {
      tree@ext <- tree@exc <- vector()
      tree@ultr <- FALSE
      tree@age <- numeric()
    }
    tree@pd <- sum(sapply(tree@ndlst, function(n) n[['spn']]))
  } else {
    tree@age <- tree@pd <- numeric()
    tree@ext <- tree@ext <- vector()
    tree@ultr <- logical()
  }
  tree@ply <- any(sapply(tree@ndlst, function(n) length(n[['ptid']]) > 2))
  initialize(tree)
}