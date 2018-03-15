#' @name calcNdBlnc
#' @title Calculate the balance of a node
#' @description Returns the balance of a node.
#' @details Balance is calculated as the absolute difference between the number of descendents 
#' of the two bifurcating edges of a node and the expected value for a balanced tree.
#' \code{NA} is returned if the node is polytomous or a tip.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{calcNdsBlnc}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' calcNdBlnc(tree, id=tree['root'])  # root balance
calcNdBlnc <- function(tree, id) {
  ntot <- length(getNdKids(tree, id))
  ptids <- tree@ndlst[[id]][['ptid']]
  if(length(ptids) > 2) {
    return(NA)
  }
  ptid <- ptids[1]
  nprt <- length(getNdKids(tree, ptid))
  if(nprt == 0) {
    nprt <- 1
  }
  abs((ntot/2) - nprt)
}

#' @name calcNdsBlnc
#' @title Calculate the balances of all nodes
#' @description Returns the absolute differences in number of descendants for bifurcating 
#' branches of every node
#' @details Runs \code{calcNdBlnc()} across all node IDs. \code{NA} is returned if the
#' node is polytomous. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids node ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{calcNdBlnc}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' calcNdsBlnc(tree, ids=tree['nds'])
calcNdsBlnc <- function(tree, ids, parallel=FALSE, progress="none") {
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  plyr::mdply(.data=l_data, .fun=calcNdBlnc, tree=tree, 
              .parallel=parallel, .progress=progress)[ ,2]
}

#' @name calcDstTrp
#' @title Calculate the triplet distance between two trees
#' @description Returns the triplet distance between two trees.
#' @details The triplet distance is calculated as the sum of different outgroups among
#' every triplet of tips between the two trees. Normalisation is performed by dividing
#' the resulting number by the total number of triplets shared between the two trees.
#' The triplet distance is calculated only for shared tips between the two trees. Parallelizable.
#' @param tree_1 \code{TreeMan} object
#' @param tree_2 \code{TreeMan} object
#' @param nrmlsd Boolean, should returned value be between 0 and 1? Default, FALSE.
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @references
#' Critchlow DE, Pearl DK, Qian C. (1996) The Triples Distance for rooted bifurcating phylogenetic trees.
#' Systematic Biology, 45, 323-34.
#' @seealso
#' \code{\link{calcDstBLD}}, \code{\link{calcDstRF}} 
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#' library(treeman)
#' tree_1 <- randTree(10)
#' tree_2 <- randTree(10)
#' calcDstTrp(tree_1, tree_2)
calcDstTrp <- function(tree_1, tree_2, nrmlsd=FALSE,
                       parallel=FALSE, progress="none") {
  .count <- function(i) {
    o1 <- getOtgrp(tree_1, cmbs[ ,i])
    o2 <- getOtgrp(tree_2, cmbs[ ,i])
    if (length(o1) != length(o2) || o1 != o2) {
      cntr <- 1
    } else {
      cntr <- 0
    }
    cntr
  }
  shrd <- tree_1@tips[tree_1@tips %in% tree_2@tips]
  cmbs <- combn(shrd, 3)
  l_data <- data.frame(i=1:ncol(cmbs), stringsAsFactors=FALSE)
  res <- plyr::mdply(.data=l_data, .count, .parallel=parallel, .progress=progress)
  cntr <- sum(res[ ,2])
  if (nrmlsd) {
    cntr <- cntr/ncol(cmbs)
  }
  cntr
}

#' @name calcOvrlp
#' @title Calculate phylogenetic overlap
#' @description Returns the sum of branch lengths represented by ids_1 and ids_2 for a tree.
#' @details Use this to calculate the sum of branch lengths that are represented between two
#' communities. This measure is also known as the unique fraction. It can be used to measure
#' concepts of phylogenetic turnover. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids_1 tip ids of community 1
#' @param ids_2 tip ids of community 2
#' @param nrmlsd Boolean, should returned value be between 0 and 1? Default, FALSE.
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @references
#' Lozupone, C., & Knight, R. (2005). UniFrac: a new phylogenetic method for comparing
#' microbial communities. Applied and Environmental Microbiology, 71(12), 8228-35.
#' @seealso
#' \code{\link{calcPhyDv}}
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' ids_1 <- sample(tree['tips'], 5)
#' ids_2 <- sample(tree['tips'], 5)
#' calcOvrlp(tree, ids_1, ids_2)
calcOvrlp <- function(tree, ids_1, ids_2, nrmlsd=FALSE,
                      parallel=FALSE, progress="none") {
  if(progress != "none") {
    cat("Part 1/3 ....\n")
  }
  spns <- getNdsSlt(tree, slt_nm='spn', ids=tree@all,
                    parallel=parallel, progress=progress)
  names(spns) <- tree@all
  if(progress != "none") {
    cat("Part 2/3 ....\n")
  }
  ids_1 <- c(unique(unlist(getNdsPrids(tree, ids=ids_1,
                                       parallel=parallel,
                                       progress=progress))), ids_1)
  if(progress != "none") {
    cat("Part 3/3 ....\n")
  }
  ids_2 <- c(unique(unlist(getNdsPrids(tree, ids=ids_2,
                                       parallel=parallel,
                                       progress=progress))), ids_2)
  ovrlp <- sum(spns[ids_2[ids_2 %in% ids_1]])
  if(nrmlsd) {
    ovrlp <- ovrlp/tree@pd
  }
  ovrlp
}

#' @name calcDstBLD
#' @title Calculate the BLD between two trees
#' @description Returns the branch length distance between two trees.
#' @details BLD is the Robinson-Foulds distance weighted by branch length. Instead of summing
#' the differences in partitions between the two trees, the metric takes the square root
#' of the squared difference in branch lengths. Parallelizable.
#' @param tree_1 \code{TreeMan} object
#' @param tree_2 \code{TreeMan} object
#' @param nrmlsd Boolean, should returned value be between 0 and 1? Default, FALSE.
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @references
#' Kuhner, M. K. and Felsenstein, J. (1994) Simulation comparison of phylogeny
#' algorithms under equal and unequal evolutionary rates. Molecular Biology and
#' Evolution, 11, 459-468.
#' @seealso
#' \code{\link{calcDstTrp}}, \code{\link{calcDstRF}} 
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#' library(treeman)
#' tree_1 <- randTree(10)
#' tree_2 <- randTree(10)
#' calcDstBLD(tree_1, tree_2)

calcDstBLD <- function(tree_1, tree_2, nrmlsd=FALSE,
                       parallel=FALSE, progress="none") {
  n1 <- tree_1@nds[!tree_1@nds == tree_1@root]
  n2 <- tree_2@nds[!tree_2@nds == tree_2@root]
  if(progress != "none") {
    cat("Part 1/2 ....\n")
  }
  c1 <- getNdsKids(tree_1, n1, parallel=parallel, progress=progress)
  if(progress != "none") {
    cat("Part 2/2 ....\n")
  }
  c2 <- getNdsKids(tree_2, n2, parallel=parallel, progress=progress)
  s1 <- getNdsSlt(tree_1, slt_nm="spn", ids=n1)
  s2 <- getNdsSlt(tree_2, slt_nm="spn", ids=n2)
  d1 <- s2[match(c1, c2)]
  d1[which(is.na(d1))] <- 0
  d1 <- s1 - d1
  d2 <- s1[match(c2, c1)]
  d2[which(is.na(d2))] <- 0
  d2 <- s2 - d2
  d <- sqrt(sum(c(d1^2, d2^2)))
  if(nrmlsd) {
    max_d <- sqrt(sum(c(s1^2, s2^2)))
    d <- d/max_d
  }
  d
}

#' @name calcDstRF
#' @title Calculate the Robinson-Foulds distance between two trees
#' @description Returns the Robinson-Foulds distance between two trees.
#' @details RF distance is calculated as the sum of partitions in one tree that are
#' not shared by the other. The maximum number of split differences is the total number
#' of nodes in both trees (excluding the roots). Trees are assumed to be bifurcating,
#' this is not tested. The metric is calculated as if trees are unrooted. Parallelizable.
#' @param tree_1 \code{TreeMan} object
#' @param tree_2 \code{TreeMan} object
#' @param nrmlsd Boolean, should returned value be between 0 and 1? Default, FALSE.
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @references
#' Robinson, D. R.; Foulds, L. R. (1981). "Comparison of phylogenetic trees".
#' Mathematical Biosciences 53: 131-147.
#' @seealso
#' \code{\link{calcDstBLD}}, \code{\link{calcDstTrp}} 
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#' library(treeman)
#' tree_1 <- randTree(10)
#' tree_2 <- randTree(10)
#' calcDstRF(tree_1, tree_2)
calcDstRF <- function(tree_1, tree_2, nrmlsd=FALSE,
                      parallel=FALSE, progress="none") {
  # get unrooted bipartitions
  rp_1 <- getNdSlt(tree=tree_1, id=tree_1@root,
                   slt_nm='ptid')
  rp_2 <- getNdSlt(tree=tree_2, id=tree_2@root,
                   slt_nm='ptid')
  rp_1 <- rp_1[!rp_1 %in% tree_1@tips]
  rp_2 <- rp_2[!rp_2 %in% tree_2@tips]
  ignr_1 <- c(rp_1[1], tree_1@root)
  ignr_2 <- c(rp_2[1], tree_2@root)
  n1 <- tree_1@nds[!tree_1@nds %in% ignr_1]
  n2 <- tree_2@nds[!tree_2@nds %in% ignr_2]
  if(progress != "none") {
    cat("Part 1/2 ....\n")
  }
  c1 <- getNdsKids(tree_1, n1, parallel=parallel,
                   progress=progress)
  c1 <- lapply(c1, sort)
  if(progress != "none") {
    cat("Part 2/2 ....\n")
  }
  c2 <- getNdsKids(tree_2, n2, parallel=parallel,
                   progress=progress)
  c2 <- lapply(c2, sort)
  d <- sum(!c1 %in% c2) + sum(!c2 %in% c1)
  if(nrmlsd) {
    max_d <- (length(n1) + length(n2))
    d <- d/max_d
  }
  d
}

#' @name calcPhyDv
#' @title Calculate phylogenetic diversity
#' @description Returns the phylogenetic diversity of a tree for the tips specified.
#' @details Faith's phylogenetic diversity is calculated as the sum of all connected
#' branches for specified tips in a tree. It can be used to investigate how biodviersity
#' as measured by the phylogeny changes. Parallelizable.
#' The function uses \code{getCnntdNds()}.
#' @param tree \code{TreeMan} object
#' @param tids tip ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @references
#' Faith, D. (1992). Conservation evaluation and phylogenetic diversity.
#'  Biological Conservation, 61, 1-10.
#' @seealso
#' \code{\link{calcFrPrp}}, \code{\link{calcOvrlp}}, \code{\link{getCnnctdNds}},
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' calcPhyDv(tree, tree['tips'])
calcPhyDv <- function(tree, tids,
                      parallel=FALSE, progress="none") {
  prids <- getCnnctdNds(tree, tids)
  spns <- getNdsSlt(tree, slt_nm="spn", ids=prids,
                    parallel=parallel, progress=progress)
  sum(spns)
}

#' @name calcFrPrp
#' @title Calculate evolutionary distinctness
#' @description Returns the evolutationary distinctness of ids using the fair proportion metric.
#' @details The fair proportion metric calculates the evolutionary distinctness of tips
#' in a tree through summing the total amount of branch length each tip represents, where
#' each branch in the tree is evenly divided between all descendants. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param tids tip IDs
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @references
#' Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. and Baillie, J.E.M. (2007). 
#'  Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE, 2, e296.
#' @seealso
#' \code{\link{calcPhyDv}}, \code{\link{calcPrtFrPrp}},
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' calcFrPrp(tree, tree['tips'])
calcFrPrp <- function(tree, tids, progress="none") {
  .calc <- function(i) {
    id <- tree@all[i]
    spn <- getNdSlt(tree, "spn", id)
    if(id %in% tree@tips) {
      spn_shres[i, id] <<- spn
    } else {
      kids <- getNdKids(tree, id)
      spn_shre <- spn/length(kids)
      spn_shres[i, kids] <<- spn_shre
    }
  }
  spn_shres <- matrix(0, ncol=tree@ntips, nrow=tree@nall)
  colnames(spn_shres) <- tree@tips
  plyr::m_ply(.data=data.frame(i=1:tree@nall), .fun = .calc,
              .progress=progress)
  colSums(spn_shres[, tids])
}

#' @name calcPrtFrPrp
#' @title Calculate evolutionary distinctness for part of tree
#' @description Returns the evolutationary distinctness of ids using the fair proportion metric.
#' @details Extension of \code{calcFrPrp()} but with ignore argument.
#' Use \code{ignr} to ignore certain tips from calculation. For example, if any of tips
#' are extinct you may wish to ignore these.
#' @param tree \code{TreeMan} object
#' @param tids tip IDs
#' @param ignr tips to ignore in calculation
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @references
#' Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. and Baillie, J.E.M. (2007). 
#'  Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE, 2, e296.
#' @seealso
#' \code{\link{calcFrPrp}}
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' calcPrtFrPrp(tree, c('t1','t3'), ignr='t2')
calcPrtFrPrp <- function(tree, tids, ignr=NULL, progress="none") {
  .calc <- function(i) {
    id <- allnds[i]
    spn <- getNdSlt(tree, "spn", id)
    if(id %in% tips) {
      spn_shres[i, id] <<- spn
    } else {
      kids <- getNdKids(tree, id)
      kids <- kids[!kids %in% ignr]
      if(length(kids) > 0) {
        spn_shre <- spn/length(kids)
        spn_shres[i, kids] <<- spn_shre
      }
    }
  }
  tips <- tree@tips[!tree@tips %in% ignr]
  allnds <- tree@all[!tree@all %in% ignr]
  spn_shres <- matrix(0, ncol=length(tips), nrow=length(allnds))
  colnames(spn_shres) <- tips
  plyr::m_ply(.data=data.frame(i=1:length(allnds)), .fun = .calc,
              .progress=progress)
  colSums(spn_shres[, tids])
}

#' @name calcDstMtrx
#' @title Calculate the distance matrix
#' @description Returns a distance matrix for specified ids of a tree.
#' @details The distance between every id in the tree is calculated by summing the
#' lengths of the branches that connect them. This can be useful for testing the distances
#' between trees, checking for evoltuionary isolated tips etc. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids IDs of nodes/tips
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{calcDstBLD}}, \code{\link{calcDstRF}}, \code{\link{calcDstTrp}}
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#' # checking the distance between two trees
#' library(treeman)
#' tree_1 <- randTree(10)
#' tree_2 <- randTree(10)
#' dmat1 <- calcDstMtrx(tree_1, tree_1['tips'])
#' dmat2 <- calcDstMtrx(tree_2, tree_2['tips'])
#' mdl <- cor.test(x=dmat1, y=dmat2)
#' as.numeric(1 - mdl$estimate)  # 1 - Pearson's r
calcDstMtrx <- function(tree, ids, parallel=FALSE,
                        progress="none") {
  .getDist <- function(id_1, id_2) {
    if(id_1 == id_2) {
      return(0)
    }
    path <- getPath(tree, from=id_1, to=id_2)
    path_spns <- unlist(lapply(tree@ndlst[path], function(n) n[['spn']]))
    sum(path_spns)
  }
  cmbs <- expand.grid(ids, ids, stringsAsFactors=FALSE)
  colnames(cmbs) <- c('id_1', 'id_2')
  res <- plyr::mdply(.data=cmbs, .fun=.getDist,
                     .parallel=parallel, .progress=progress)
  res <- matrix(res[ ,3], ncol=length(ids))
  colnames(res) <- rownames(res) <- ids
  res
}