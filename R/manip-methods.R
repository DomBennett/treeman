# TODO: mergeTree, collapseNode, removeNode, reRoot

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
  if(tree@updtd) {
    tree@ndmtrx <- tree@ndmtrx[bool, bool]
    tree <- updateTree(tree)
  }
  tree
}

# TODO: fix
#' #' @name addTips
#' #' @title Add tips to a tree
#' #' @description Returns a tree with tip ID(s) added
#' #' @details User must provide new tip ID(s), the ID(s) of the nodes
#' #' which will become the new tip's sister, and new branch lengths.
#' #' Optionally, user can specify the IDs for the new parent nodes.
#' #' @param tree \code{TreeMan} object
#' #' @param tids tip IDs
#' #' @param sids IDs of node that will become new tip sisters
#' #' @param spns length of new branch
#' #' @param pids parent ID (default is 'p_' + tid)
#' #' @seealso
#' #' \code{\link{rmTips}}, 
#' #' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' #' @export
#' #' @examples
#' #' library(treeman)
#' #' tree <- randTree(10)
#' #' tree <- addTip(tree, tids='t11', sids='t1', spns=1)
#' #' tree <- updateTree(tree)
#' #' summary(tree)
#' addTip <- function(tree, tids, sids, pids=paste0("p_", tids)) {
#'   # internals
#'   .addTip <- function(i) {
#'     # terminology
#'     # snd, sid -- old sister node and id
#'     # tnd, tid -- new tip node and id
#'     # pnd, pid -- new parent node and id
#'     # gpnd, gpid -- grand parent (prid of old sister)
#'     # unpack
#'     tid <- tids[i]
#'     sids <- sids[i]
#'     pid <- pids[i]
#'     spn <- spns[i]
#'     # init new nodes
#'     tnd <- list('id'=tid)
#'     snd <- ndlst[[sid]]
#'     gpid <- snd[['prid']][[1]]
#'     gpnd <- ndlst[[gpid]]
#'     pnd <- list('id'=pid, 'kids'=sid)
#'     # update spans
#'     tnd[['spn']] <- start - end
#'     pnd[['spn']] <- snd[['spn']] - (start - age)
#'     snd[['spn']] <- start - age
#'     # update ptid
#'     gpnd[['ptid']] <- gpnd[['ptid']][!gpnd[['ptid']] %in% snd[['id']]]
#'     gpnd[['ptid']] <- c(gpnd[['ptid']], pnd[['id']])
#'     pnd[['ptid']] <- c(tid, sid)
#'     # set prid
#'     tnd[['prid']] <- snd[['prid']]
#'     pnd[['prid']] <- snd[['prid']]
#'     # add to ndlst
#'     ndlst[[tid]] <<- tnd
#'     ndlst[[pid]] <<- pnd
#'     ndlst[[sid]] <<- snd
#'     ndlst[[gpid]] <<- gpnd
#'   }
#'   ndlst <- tree@ndlst
#'   sapply(1:length(ids))
#'   tree@ndlst <- ndlst
#'   tree
#' }
#' 
#' #' @name pinTips
#' #' @title Pin tips to a tree
#' #' @description Returns a tree with new tips added based on given lineages and time points
#' #' @details User must provide a vector of new tip IDs, a list of the ranked lineages
#' #' for these IDs (in ascending order) and a vector of end time points for each new ID
#' #' (0s for extant tips). The function expects the given tree to be taxonomically informed;
#' #' the \code{txnym} slot for every node should have a taxonomic label. The function takes
#' #' the lineage and tries to randomly add the new tip at the lowest point in the taxonomic rank
#' #' before the end time point. Parallelizable.
#' #' @param tree \code{TreeMan} object
#' #' @param tids new tip ids
#' #' @param lngs list of vectors of the lineages of each tid
#' #' @param ends end time points for each tid
#' #' @param ... \code{plyr} arguments
#' #' @seealso
#' #' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' #' @export
#' #' @examples
#' #' # see https://github.com/DomBennett/treeman/wiki/Pinning-tips for a detailed example
#' pinTips <- function(tree, tids, lngs, ends, ...) {
#'   .pin <- function(i) {
#'     # unpack
#'     tid <- tids[i]
#'     lng <- lngs[[i]]
#'     end <- ends[i]
#'     for(j in length(lng):1) {
#'       spns <- names(txnyms)[which(txnyms %in% lng[j])]
#'       if(length(spns) == 0) {
#'         next
#'       }
#'       spns <- c(spns, unlist(sapply(spns, function(n) tree@ndlst[[n]][['ptid']])))
#'       spns <- spns[spns != tree@root]
#'       rngs <- getSpnsAge(tree, ids=spns)
#'       bool <- rngs[ ,'start'] > end
#'       if(any(bool)) {
#'         rngs <- rngs[bool, ]
#'         rngs[rngs[ ,'end'] <= end, "end"] <- end
#'         # pinning is based on branch length
#'         prbs <- rngs$start - rngs$end
#'         e <- as.vector(sample(rngs$spn, prob=prbs, size=1))
#'         e_i <- which(rngs$spn == e)
#'         start <- runif(min=rngs$end[e_i], max=rngs$start[e_i], n=1)
#'         if(j != length(lng)) {
#'           tip_txnym <- lng[j+1]
#'         } else {
#'           tip_txnym <- lng[j]
#'         }
#'         pid <- paste0('p_', tid, sep='')
#'         tree <- addTip(tree, tid=tid, sid=e, start=start, end=end,
#'                         pid=pid)
#'         tree@ndlst[[tid]][['txnym']] <- tip_txnym
#'         tree@ndlst[[pid]][['txnym']] <- lng[j]
#'         # add to txnyms list
#'         txnyms[[tid]] <<- tip_txnym
#'         # push out
#'         tree <<- tree
#'         break
#'       }
#'     }
#'   }
#'   .getTxnyms <- function(txnym, ...) {
#'     txnym
#'   }
#'   txnyms <- plyr::mlply(tree@ndlst, .fun=.getTxnyms)
#'   txnyms <- txnyms[1:length(txnyms)]
#'   names(txnyms) <- names(tree@ndlst)
#'   plyr::m_ply(1:length(tids), .pin)
#'   tree
#' }
