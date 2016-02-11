# TODO: addTips, removeTip, mergeTree, collapseNode, removeNode
# TODO: add doc for adding and removing tips

rmTip <- function(...) {
  cat('Yep... someone needs to create this function. Sorry!\n')
}

addTip <- function(tree, id, sister, start, end,
                   parent_id=paste0("p_", id),
                   tip_txnym=NULL, parent_txnym=NULL) {
  tip <- list('id'=id)
  if(!is.null(tip_txnym)) {
    tip[['txnym']] <- tip_txnym
  }
  node <- list('id'=parent_id)
  if(!is.null(parent_txnym)) {
    node[['txnym']] <- parent_txnym
  }
  tip[['span']] <- start - end
  age <- getNodeAge(tree, sister)
  new_sister <- sister <- tree@nodelist[[sister]]
  new_parent <- tree@nodelist[[sister[['prid']][1]]]
  new_parent[['ptid']] <- new_parent[['ptid']][!new_parent[['ptid']] %in% sister[['id']]]
  new_parent[['ptid']] <- c(new_parent[['ptid']], node[['id']])
  new_sister[['span']] <- start - age
  new_sister[['prid']] <- c(node[['id']], sister[['prid']])
  node[['span']] <- sister[['span']] - new_sister[['span']]
  node[['pd']] <- new_sister[['span']] + tip[['span']]
  node[['prdst']] <- sister[['prdst']] - new_sister[['span']]
  node[['prid']] <- sister[['prid']]
  node[['ptid']] <- node[['kids']] <- c(tip[['id']], sister[['id']])
  tip[['pd']] <- 0
  tip[['prdst']] <- node[['prdst']] + tip[['span']]
  tip[['prid']] <- new_sister[['prid']]
  ndlst <- tree@nodelist
  tree@nodelist[[tip[['id']]]] <- tip
  tree@nodelist[[node[['id']]]] <- node
  tree@nodelist[[new_sister[['id']]]] <- new_sister
  tree@nodelist[[new_parent[['id']]]] <- new_parent
  tree@nodelist <- treeman:::.updateTip(tree@nodelist, tid=id, rid=tree@root)
  treeman:::.updateSlots(tree)
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
        prbs <- rngs$start - rngs$end
        e <- as.vector(sample(rngs$span, prob=prbs, size=1))
        e_i <- which(rngs$span == e)
        start <- runif(min=rngs$end[e_i], max=rngs$start[e_i], n=1)
        if(j != length(lng)) {
          tip_txnym <- lng[j+1]
        } else {
          tip_txnym <- lng[j]
        }
        tree <<- addTip(tree, id=tid, sister=e, start=start, end=end,
                       tip_txnym=tip_txnym, parent_txnym=lng[j])
        # add to txnyms list
        txnyms[[tid]] <<- tip_txnym
        break
      }
    }
  }
  .getTxnyms <- function(txnym, ptid, ...) {
    if(!exists('ptid')) {
      # if tip node, only take first element
      res <- txnym[1]
    } else {
      res <- txnym
    }
    res
  }
  txnyms <- mlply(tree@nodelist, .fun=.getTxnyms)
  txnyms <- txnyms[1:length(txnyms)]
  names(txnyms) <- tree@all
  m_ply(1:length(tids), .pin)
  tree
}
