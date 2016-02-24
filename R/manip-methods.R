# TODO: addTips, removeTip, mergeTree, collapseNode, removeNode
# TODO: add doc for adding and removing tips

rmTip <- function(...) {
  cat('Yep... someone needs to create this function. Sorry!\n')
}

addTip <- function(tree, tid, sid, start, end,
                   pid=paste0("p_", tid)) {
  # terminology
  # snd, sid -- old sister node and id
  # tnd, tid -- new tip node and id
  # pnd, pid -- new parent node and id
  # gpnd, gpid -- grand parent (prid of old sister)
  updater <- function(nd) {
    # node operation
    i <- which(nd[['prid']] == gpid) - 1
    mxi <- length(nd[['prid']])
    nd[['prid']] <- c(nd[['prid']][0:i], pid, 
                      nd[['prid']][(i+1):mxi])
    nd
  }
  # unpack
  ndlst <- tree@nodelist
  # get key data from tree
  nids <- getNodePtid(tree, sid)
  age <- getNodeAge(tree, sid)
  # init new nodes
  tnd <- list('id'=tid)
  snd <- ndlst[[sid]]
  gpid <- snd[['prid']][[1]]
  gpnd <- ndlst[[gpid]]
  pnd <- list('id'=pid, 'kids'=sid)
  # update spans
  tnd[['span']] <- start - end
  pnd[['span']] <- snd[['span']] - (start - age)
  snd[['span']] <- start - age
  # update ptid
  gpnd[['ptid']] <- gpnd[['ptid']][!gpnd[['ptid']] %in% snd[['id']]]
  gpnd[['ptid']] <- c(gpnd[['ptid']], pnd[['id']])
  pnd[['ptid']] <- c(tid, sid)
  # set prid
  tnd[['prid']] <- snd[['prid']]
  pnd[['prid']] <- snd[['prid']]
  # set prdst
  pnd[['prdst']] <- snd[['prdst']] - snd[['span']]
  tnd[['prdst']] <- pnd[['prdst']] + tnd[['span']]
  # set kids
  if(is.null(snd[['kids']])) {
    pnd[['kids']] <- sid
  } else {
    pnd[['kids']] <- snd[['kids']]
  }
  # set pd
  pnd[['pd']] <- snd[['pd']]
  tnd[['pd']] <- 0
  # add to ndlst
  ndlst[[tid]] <- tnd
  ndlst[[pid]] <- pnd
  ndlst[[sid]] <- snd
  ndlst[[gpid]] <- gpnd
  # update upstream prids from sid onwards
  ndlst <- .updateNodesSlot(ndlst, c(tid, nids), updater)
  # update downstream using updateTip
  tree@nodelist <- .updateTip(ndlst, tid=tid, rid=tree@root)
  .updateTreeSlots(tree)
}

pinTips <- function(tree, tids, lngs, ends) {
  .pin <- function(i) {
    # unpack
    tid <- tids[i]
    lng <- lngs[[i]]
    end <- ends[i]
    for(j in length(lng):1) {
      spns <- names(txnyms)[which(txnyms %in% lng[j])]
      if(length(spns) == 0) {
        next
      }
      spns <- c(spns, unlist(sapply(spns, function(n) tree@nodelist[[n]][['ptid']])))
      spns <- spns[spns != tree@root]
      rngs <- getSpansAge(tree, ids=spns)
      bool <- rngs[ ,'start'] > end
      if(any(bool)) {
        rngs <- rngs[bool, ]
        rngs[rngs[ ,'end'] <= end, "end"] <- end
        # pinning is based on branch length
        prbs <- rngs$start - rngs$end
        e <- as.vector(sample(rngs$span, prob=prbs, size=1))
        e_i <- which(rngs$span == e)
        start <- runif(min=rngs$end[e_i], max=rngs$start[e_i], n=1)
        if(j != length(lng)) {
          tip_txnym <- lng[j+1]
        } else {
          tip_txnym <- lng[j]
        }
        pid <- paste0('p_', tid, sep='')
        tree <<- addTip(tree, tid=tid, sid=e, start=start, end=end,
                        pid=pid)
        tree@nodelist[[tid]][['taxonym']] <- tip_txnym
        tree@nodelist[[pid]][['taxonym']] <- lng[j]
        # add to txnyms list
        txnyms[[tid]] <<- tip_txnym
        break
      }
    }
  }
  .getTxnyms <- function(txnym, ...) {
    txnym
  }
  txnyms <- mlply(tree@nodelist, .fun=.getTxnyms)
  txnyms <- txnyms[1:length(txnyms)]
  names(txnyms) <- names(tree@nodelist)
  m_ply(1:length(tids), .pin)
  tree
}
