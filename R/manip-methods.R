# TODO: mergeTree, collapseNode, removeNode

#' @name rmTips
#' @title Remove tips from a tree
#' @description Returns a tree with a tip ID(s) for removal
#' @details Removes tips in a tree. Set drp_intrnl to FALSE to convert
#' internal nodes into new tips.
#' @param tree \code{TreeMan} object
#' @param tids tip IDs
#' @param drp_intrnl Boolean, drop internal branches, default FALSE
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' tree <- rmTips(tree, 't1')
#' summary(tree)
# @seealso
# \code{\link{addTip}}, 
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
#' Optionally, user can specify the IDs for the new parental interal nodes.
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
#' @param lngs list of vectors of the lineages of each tid
#' @param end_ages end time points for each tid
#' @param tree_age age of tree
#' @seealso
#' \code{\link{addTip}}, \code{\link{rmTips}},
#' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#' # see https://github.com/DomBennett/treeman/wiki/Pinning-tips for a detailed example
pinTips <- function(tree, tids, lngs, end_ages, tree_age) {
  .pin <- function(i) {
    # unpack
    tid <- tids[i]
    lng <- lngs[[i]]
    end <- end_ages[i]
    for(j in length(lng):1) {
      spns <- names(txnyms)[which(txnyms %in% lng[j])]
      if(length(spns) == 0) {
        next
      }
      spns <- unique(c(spns, unlist(getNdsPtids(tree, spns))))
      spns <- spns[spns != tree@root]
      rngs <- getSpnsAge(tree, ids=spns, tree_age=tree_age)
      bool <- rngs[ ,'start'] > end
      if(any(bool)) {
        rngs <- rngs[bool, ]
        rngs[rngs[ ,'end'] <= end, "end"] <- end
        # pinning is based on branch length
        prbs <- rngs$start - rngs$end
        e <- as.vector(sample(rngs$spn, prob=prbs, size=1))
        e_i <- which(rngs$spn == e)
        start <- runif(min=rngs$end[e_i], max=rngs$start[e_i], n=1)
        tip_txnym <- lng[length(lng)]
        if(j == length(lng)) {
          pid_txnym <- lng[length(lng)]
        } else {
          pid_txnym <- lng[j:length(lng)]
        }
        pid <- paste0('p_', tid, sep='')
        tree <- addTip(tree, tid=tid, sid=e, strt_age=start, end_age=end,
                       pid=pid, tree_age=tree_age)
        tree@ndlst[[tid]][['txnym']] <- tip_txnym
        tree@ndlst[[pid]][['txnym']] <- pid_txnym
        # add to txnyms list
        txnyms[[tid]] <<- tip_txnym
        txnyms[[pid]] <<- pid_txnym
        # push out
        tree <<- tree
        break
      }
    }
  }
  .testLngs <- function(lng) {
    for(l in lng) {
      if(grepl('[^a-zA-Z_0-9]', l)) {
        stop(paste0('Unsuitable characters in [', l, ']'))
      }
    }
    NULL
  }
  sapply(lngs, .testLngs)
  .getTxnyms <- function(txnym, ...) {
    txnym
  }
  if(!tree@wtxnyms) {
    stop('tree has no txnyms')
  }
  if(any(end_ages < 0)) {
    warning('One or more end ages are less than zero, this will change the age of tree.')
  }
  txnyms <- plyr::mlply(tree@ndlst, .fun=.getTxnyms)
  txnyms <- txnyms[1:length(txnyms)]
  names(txnyms) <- names(tree@ndlst)
  plyr::m_ply(1:length(tids), .pin)
  tree
}
