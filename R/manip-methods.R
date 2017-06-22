
#' @name ultrTree
#' @title Make tree ultrametric
#' @description Returns a tree with all tips ending at time 0
#' @details Re-calculates the branch lengths in the tree so that all
#' tips are brought to the same time point: all species are extant.
#' @param tree \code{TreeMan} object
#' @seealso
#' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' (getDcsd(tree))  # list all extinct tips
#' tree <- ultrTree(tree)
#' (getDcsd(tree))  # list all extinct tips
ultrTree <- function(tree) {
  # bring tips to maximum possible length
  tip_dpths <- sapply(getNdsPrids(tree, tree@tips), length)
  tip_spns <- (tree['ntips'] - tip_dpths)
  nd_spns <- rep(1, tree@nnds)
  names(nd_spns) <- tree@nds
  spns <- c(tip_spns, nd_spns)
  tree <- setNdsSpn(tree, ids=names(spns), vals=spns)
  updateSlts(tree)
}

#' @name rmNodes
#' @title Remove nodes from a tree
#' @description Returns a tree with a node ID(s) removed
#' @details Removes nodes in a tree. Joins the nodes following to
#' the nodes preceding the node to be removed. Creates polytomies.
#' Warning: do not use this function to remove tip nodes, this create a
#' corrupted tree.
#' @param tree \code{TreeMan} object
#' @param nids internal node IDs
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{addTip}}, \code{\link{rmTips}},
#'  \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' tree <- rmNodes(tree, 'n3')
#' summary(tree)  # tree is now polytmous
rmNodes <- function(tree, nids, progress='none') {
  .rmNode <- function(nid) {
    ptids <- ndlst[[nid]][['ptid']]
    prid <- ndlst[[nid]][['prid']]
    if(tree@wspn) {
      for(ptid in ptids) {
        ndlst[[ptid]][['spn']] <-
          ndlst[[ptid]][['spn']] +
          ndlst[[nid]][['spn']]
      }
    }
    for(ptid in ptids) {
      ndlst[[ptid]][['prid']] <- prid
    }
    new_ptids <- ndlst[[prid]][['ptid']]
    new_ptids <- new_ptids[new_ptids != nid]
    new_ptids <- c(new_ptids, ptids)
    ndlst[[prid]][['ptid']] <- new_ptids
    ndlst <<- ndlst[names(ndlst) != nid]
  }
  if(tree@root %in% nids) {
    stop('Cannot remove root.')
  }
  ndlst <- tree@ndlst
  plyr::m_ply(.data=nids, .fun=.rmNode, .progress=progress)
  bool <- tree@all %in% names(ndlst)
  tree@ndlst <- ndlst
  tree <- pstMnp(tree)
  tree <- updateSlts(tree)
  if(!is.null(tree@ndmtrx)) {
    tree@ndmtrx <- bigmemory::as.big.matrix(tree@ndmtrx[bool, bool])
  }
  tree
}

#' @name rmTips
#' @title Remove tips from a tree
#' @description Returns a tree with a tip ID(s) removed
#' @details Removes tips in a tree. Set drp_intrnl to FALSE to convert
#' internal nodes into new tips. Warning: do not use this function to remove
#' internal nodes, this create a corrupted tree.
#' @param tree \code{TreeMan} object
#' @param tids tip IDs
#' @param drp_intrnl Boolean, drop internal branches, default FALSE
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{addTip}}, \code{\link{rmNodes}},
#'  \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' tree <- rmTips(tree, 't1')
#' summary(tree)
#' # running the function using an internal
#' # node will create a corrupted tree
#' tree <- rmTips(tree, 'n3')
#' # run summary() to make sure a change has
#' # not created a corruption
#' #summary(tree)
rmTips <- function(tree, tids, drp_intrnl=TRUE, progress="none") {
  # internal
  .rmTip <- function(tid) {
    # get sister IDs
    sids <- .getNdSstrFrmLst(ndlst, tid)
    # get prid
    prid <- ndlst[[tid]][['prid']][[1]]
    # remove tid
    ndlst <<- ndlst[names(ndlst) != tid]
    ndlst[[prid]][['ptid']] <<-
      ndlst[[prid]][['ptid']][ndlst[[prid]][['ptid']] != tid]
    # remove prnd if specified and not polytomous
    if(drp_intrnl & length(sids) == 1) {
      ptid <- ndlst[[prid]][['ptid']][[1]]
      ndlst[[ptid]][['spn']] <<- ndlst[[prid]][['spn']] +
        ndlst[[ptid]][['spn']]
      if(prid != rid) {
        gprid <- ndlst[[prid]][['prid']][[1]]
        ndlst[[ptid]][['prid']] <<- gprid
        g_ptids <- ndlst[[gprid]][['ptid']]
        g_ptids <- g_ptids[g_ptids != prid]
        ndlst[[gprid]][['ptid']] <<- c(g_ptids, ptid)
      } else {
        # if prid to be dropped is root, set root to ptid
        rid <<- ptid
        ndlst[[ptid]][['prid']] <<- ptid
      }
      ndlst <<- ndlst[names(ndlst) != prid]
    }
  }
  ndlst <- tree@ndlst
  rid <- tree@root
  plyr::m_ply(.data=tids, .fun=.rmTip, .progress=progress)
  bool <- tree@all %in% names(ndlst)
  tree@ndlst <- ndlst
  tree@root <- rid
  tree <- pstMnp(tree)
  tree <- updateSlts(tree)
  if(!is.null(tree@ndmtrx)) {
    tree@ndmtrx <- bigmemory::as.big.matrix(tree@ndmtrx[bool, bool])
  }
  tree
}

#' @name addTip
#' @title Add tip to a tree
#' @description Returns a tree with a new tip ID added
#' @details User must provide new tip ID, the ID of the node
#' which will become the new tip's sister, and new branch lengths.
#' The tip ID must only contain letters numbers and underscores.
#' Optionally, user can specify the IDs for the new parental internal nodes.
#' Ensure that the \code{strt_age} is greater than the \code{end_age}, and that
#' the \code{strt_age} falls within the age span of the sister ID. Otherwise, negative
#' spns may be produced leading to an error.
#' Note, returned tree will not have a node matrix.
#' Note, providing negative end ages will increase the age of the tree.
#' @param tree \code{TreeMan} object
#' @param tid tip ID
#' @param sid ID of node that will become new tip sisters
#' @param strt_age timepoint at which new tips first appear in the tree
#' @param end_age timepoint at which new tips end appear in the tree, default 0.
#' @param tree_age age of tree
#' @param pid parent ID (default is 'p_' + tid)
#' @seealso
#' \code{\link{rmTips}},
#' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' tree_age <- getAge(tree)
#' possible_ages <- getSpnAge(tree, 't1', tree_age)
#' start_age <- runif(1, possible_ages[['end']], possible_ages[['start']])
#' end_age <- possible_ages[['end']]
#' tree <- addTip(tree, tid='t11', sid='t1', strt_age=start_age,
#' end_age=end_age, tree_age=tree_age)
#' summary(tree)
addTip <- function(tree, tid, sid, strt_age=NULL,
                   end_age=0, tree_age=NULL,
                   pid=paste0("p_", tid)) {
  if(!is.numeric(strt_age) & tree@wspn) {
    stop('Valid strt_age not provided')
  }
  ndlst <- tree@ndlst
  # terminology
  # snd, sid -- old sister node and id
  # tnd, tid -- new tip node and id
  # pnd, pid -- new parent node and id
  # gpnd, gpid -- grand parent (prid of old sister)
  if(grepl('[^a-zA-Z_0-9]', tid)) {
    stop(paste0('Unsuitable characters in tid [', tid, ']'))
  }
  # init new nodes
  tnd <- list('id'=tid, 'prid'=pid, 'ptid'=character(), 'spn'=0)
  snd <- ndlst[[sid]]
  gpid <- snd[['prid']][[1]]
  gpnd <- ndlst[[gpid]]
  pnd <- list('id'=pid, 'spn'=0)
  # update ptid
  gpnd[['ptid']] <- gpnd[['ptid']][!gpnd[['ptid']] %in% snd[['id']]]
  gpnd[['ptid']] <- c(gpnd[['ptid']], pnd[['id']])
  pnd[['ptid']] <- c(tid, sid)
  # set prid
  pnd[['prid']] <- snd[['prid']]
  snd[['prid']] <- pid
  # add to ndlst
  ndlst[[tid]] <- tnd
  ndlst[[pid]] <- pnd
  ndlst[[sid]] <- snd
  ndlst[[gpid]] <- gpnd
  if(tree@wspn) {
    sspn <- ndlst[[sid]][['spn']]
    prid <- ndlst[[tid]][['prid']]
    gprid <- ndlst[[prid]][['prid']]
    # calc
    gp_age <- getNdAge(tree, gprid, tree_age)
    tspn <- strt_age - end_age
    pspn <- gp_age - strt_age
    new_spsn <- abs(sspn - pspn)
    if(tspn < 0 | pspn < 0 | new_spsn < 0) {
      stop('Invalid ages given: negative spns')
    }
    ndlst[[tid]][['spn']] <- tspn
    ndlst[[prid]][['spn']] <- pspn
    ndlst[[sid]][['spn']] <- new_spsn
  }
  tree@ndlst <- ndlst
  tree <- pstMnp(tree)
  tree <- rmNdmtrx(tree)
  tree <- updateSlts(tree)
  tree
}

#' @name addClade
#' @title Add clade to tree
#' @description Returns a tree with added clade
#' @details Add a \code{TreeMan} object to an existing \code{TreeMan}
#' object by specifying an ID at which to attach. If the id specified
#' is an internal node, then the original clade descending from that
#' node will be replaced. Before running, ensure no IDs are shared
#' between the \code{tree} and the \code{clade}, except for the IDs in the clade
#' of that tree that will be replaced.
#' Note, returned tree will not have a node matrix.
#' @param tree \code{TreeMan} object
#' @param id tip/node ID in tree to which the clade will be added
#' @param clade \code{TreeMan} object
#' @seealso
#' \code{\link{rmClade}}, \code{\link{getSubtree}},
#' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#' library(treeman)
#' t1 <- randTree(100)
#' # extract a clade
#' cld <- getSubtree(t1, 'n2')
#' # remove the same clade
#' t2 <- rmClade(t1, 'n2')
#' # add the clade again
#' t3 <- addClade(t2, 'n2', cld)
#' # t1 and t3 should be the same
#' # note there is no need to remove a clade before adding
#' t3 <- addClade(t1, 'n2', cld)  # same tree
addClade <- function(tree, id, clade) {
  if(!id %in% tree@tips) {
    tree <- rmClade(tree, id)
  }
  cld_ids <- names(clade@ndlst)[names(clade@ndlst) != clade@root]
  if(any(names(tree@ndlst) %in% cld_ids) &
     any(cld_ids %in% names(tree@ndlst))) {
    stop('IDs in `clade` exist in parts of `tree` not to be replaced.',
         ' Consider checking for duplicates or renaming IDs in either `tree` or `clade`')
  }
  cld_ptids <- clade@ndlst[[clade@root]][['ptid']]
  for(cld_ptid in cld_ptids) {
    clade@ndlst[[cld_ptid]][['prid']] <- id
  }
  cld_ndlst <- clade@ndlst[names(clade@ndlst) != clade@root]
  tree@ndlst[[id]][['ptid']] <- cld_ptids
  tree@ndlst <- c(tree@ndlst, cld_ndlst)
  tree <- pstMnp(tree)
  tree <- rmNdmtrx(tree)
  updateSlts(tree)
}

#' @name rmClade
#' @title Remove a clade from a tree
#' @description Returns a tree with a clade removed
#' @details Inverse function of \code{getSubtree()}. Takes a tree
#' and removes a clade based on an internal node specified. Node
#' is specified with \code{id}, all descending nodes and tips are removed.
#' The resulting tree will replace the missing clade with a tip of \code{id}.
#' @param tree \code{TreeMan} object
#' @param id node ID parent of clade to be removed
#' @seealso
#' \code{\link{addClade}}, \code{\link{getSubtree}}, \code{\link{rmTips}}
#' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#' library(treeman)
#' t1 <- randTree(100)
#' # remove a clade
#' t2 <- rmClade(t1, 'n2')
#' summary(t1)
#' summary(t2)
rmClade <- function(tree, id) {
  ptids <- getNdPtids(tree, id)
  bool <- !tree@all %in% ptids
  tree@ndlst <- tree@ndlst[!names(tree@ndlst) %in% ptids]
  tree@ndlst[[id]][['ptid']] <- character()
  tree <- pstMnp(tree)
  if(!is.null(tree@ndmtrx)) {
    tree@ndmtrx <- bigmemory::as.big.matrix(tree@ndmtrx[bool, bool])
  }
  updateSlts(tree)
}

#' @name pinTips
#' @title Pin tips to a tree
#' @description Returns a tree with new tips added based on given lineages and time points
#' @details User must provide a vector of new tip IDs, a list of the ranked lineages
#' for these IDs (in ascending order) and a vector of end time points for each new ID
#' (0s for extant tips). The function expects the given tree to be taxonomically informed;
#' the \code{txnym} slot for every node should have a taxonomic label. The function takes
#' the lineage and tries to randomly add the new tip at the lowest point in the taxonomic rank
#' before the end time point. Note, returned tree will not have a node matrix.
#' @param tree \code{TreeMan} object
#' @param tids new tip ids
#' @param lngs list of vectors of the lineages of each tid (ordered high to low rank)
#' @param end_ages end time points for each tid
#' @param tree_age age of tree
#' @seealso
#' \code{\link{addTip}}, \code{\link{rmTips}},
#' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#' # see https://github.com/DomBennett/treeman/wiki/Pinning-tips for a detailed example
pinTips <- function(tree, tids, lngs, end_ages, tree_age) {
  .getPtntls <- function(lng, end) {
    sccs <- FALSE
    # loop backwards through taxonomy
    # genus --> family --> ....
    for(i in length(lng):1) {
      pull <- txnyms %in% lng[i]
      if(sum(pull) == 0) {
        next
      }
      ptntls <- names(txnyms)[pull]
      if(length(ptntls) == 1) {
        prnt <- ptntls
      } else {
        # assumes monophylly
        prnt <- ptntls[which.max(spn_data[ptntls, 'end'])]
      }
      ptntls <- c(prnt, getNdPtids(tree, prnt))
      ptntls <- ptntls[ptntls != rid]
      pull <- spn_data[ptntls, 'start'] > end
      if(any(pull)) {
        ptntls <- ptntls[pull]
        sccs <- TRUE
        break
      }
    }
    if(sccs) {
      return(ptntls)
    }
    NULL
  }
  .getPTxnym <- function(tip_txnym, sid) {
    gp_txnym <- txnyms[[getNdSlt(tree, 'prid', sid)]]
    s_txnym <- txnyms[[sid]]
    if(s_txnym == tip_txnym) {
      pid_txnym <- tip_txnym
    } else {
      pid_txnym <- gp_txnym
    }
    pid_txnym
  }
  .pin <- function(i) {
    tid <- tids[i]
    end <- end_ages[i]
    lng <- lngs[[i]]
    ptntls <- .getPtntls(lng, end)
    if(is.null(ptntls)) {
      message(paste0('[', tid, '] could not be added'))
      return(NULL)
    }
    rngs <- spn_data[ptntls, , drop=FALSE]
    rngs[rngs[,'end'] <= end, 'end'] <- end
    # pinning is based on branch length
    # this is not a model, it just ensures
    # taxonomically matching branch lengths
    # of the tree have equal chance.
    prbs <- rngs[ ,'start'] - rngs[ ,'end']
    if(sum(prbs) == 0) {
      message(paste0('[', tid, '] could not be added'))
      return(NULL)
    }
    sid <- as.vector(sample(ptntls, prob=prbs, size=1))
    start <- runif(min=rngs[sid, 'end'], max=rngs[sid, 'start'], n=1)
    # taxnomy of tip and parent tip based grandparent
    tip_txnym <- lng[length(lng)]
    pid_txnym <- .getPTxnym(tip_txnym, sid)
    pid <- paste0('p_', tid, sep='')
    # add tip
    tree <- addTip(tree, tid=tid, sid=sid, strt_age=start,
                   end_age=end, pid=pid, tree_age=tree_age)
    tree@ndlst[[tid]][['txnym']] <- tip_txnym
    tree@ndlst[[pid]][['txnym']] <- pid_txnym
    # add to spn_data
    tid_spn <- getSpnAge(tree, tid, tree_age)
    spn_data[tid, 'start'] <<- tid_spn[ ,'start']
    spn_data[tid, 'end'] <<- tid_spn[ ,'end']
    pid_spn <- getSpnAge(tree, pid, tree_age)
    spn_data[pid, 'start'] <<- pid_spn[ ,'start']
    spn_data[pid, 'end'] <<- pid_spn[ ,'end']
    sid_spn <- getSpnAge(tree, sid, tree_age)
    spn_data[sid, 'start'] <<- sid_spn[ ,'start']
    spn_data[sid, 'end'] <<- sid_spn[ ,'end']
    # add to txnyms list
    txnyms[[tid]] <<- tip_txnym
    txnyms[[pid]] <<- pid_txnym
    # push tree out
    tree <<- tree
  }
  .testLngs <- function(lng) {
    for(l in lng) {
      if(grepl('[^a-zA-Z_0-9]', l)) {
        stop(paste0('Unsuitable characters in [', l, ']'))
      }
    }
    NULL
  }
  if(!tree@wtxnyms) {
    stop('tree has no txnyms')
  }
  if(any(end_ages < 0)) {
    warning('One or more end ages are less than zero, this will change the age of tree.')
  }
  # make sure lineages are right
  sapply(lngs, .testLngs)
  # generate taxonomy and span data
  txnyms <- getNdsSlt(tree, 'txnym', tree@all)
  txnyms <- c(txnyms, rep(NA, length(tids)*2))
  names(txnyms) <- c(names(tree@ndlst), tids, paste0('p_', tids))
  spn_data <- matrix(NA, nrow=(length(tids)*2) + tree@nall,
                     ncol=2)
  colnames(spn_data) <- c('start', 'end')
  tmp_spn_data <- getSpnsAge(tree, tree@all, tree_age)
  rownames(spn_data) <- c(tree@all, tids, paste0('p_', tids))
  spn_data[tree@all, 'start'] <- tmp_spn_data[['start']]
  spn_data[tree@all, 'end'] <- tmp_spn_data[['end']]
  rm(tmp_spn_data)
  rid <- tree@root
  # add oldest to youngest
  ordrd <- order(end_ages, decreasing=TRUE)
  plyr::m_ply(ordrd, .pin)
  tree
}
