#' @name writeTree
#' @title Write a Newick tree
#' @description Creates a Newick tree from a \code{TreeMan} object.
#' @details The \code{ndLabels} argument can be used to add a user defined node label in
#' the Newick tree. It should take only 1 argument, \code{nd}, the node represented as a list.
#' It should only return a single character value that can be added to a newick string.
#' @param tree \code{TreeMan} object
#' @param file file path
#' @param ndLabels node label function
#' @seealso
#' \code{\link{readTree}}, \code{\link{randTree}}, \url{https://en.wikipedia.org/wiki/Newick_format}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' ndLabels <- function(n) {
#' paste0(n[['id']], '_ndlabel')
#' }
#' writeTree(tree, file='example.tre', ndLabels)
#' file.remove('example.tre')
# TODO: test this with unrooted trees, adapt for TreeMen
writeTree <- function(tree, file, ndLabels=function(nd){
  return(NULL)
  }) {
  tipBytip <- function(i) {
    ids <- c(ndlst[[prid]][['kids']], prid,
             ndlst[[prid]][['prid']])
    id <<- ids[!ids %in% deja_vues][1]
    deja_vues[i] <<- id
    spn <- ndlst[[id]][['spn']]
    if(id %in% tids) {
      prids <- getNdPrids(tree, id)
      dpth <- which(prids == prid) - 1
      prid <<- ndlst[[id]][['prid']]
      tpstr <- paste0(id, ':', spn)
      if(dpth > 0) {
        brckts <- paste0(rep('(', dpth), collapse='')
        trstr <<- paste0(trstr, ',', brckts, tpstr)
      } else {
        trstr <<- paste0(trstr, ',', tpstr)
      }
    } else {
      prid <<- ndlst[[id]][['prid']]
      ndlbl <- ndLabels(ndlst[[id]])
      trstr <<- paste0(trstr, ')', ndlbl,':', spn)
    }
    NULL
  }
  # start with first tip
  # loop through tree structure adding tip by tip to string
  # unpack
  ndlst <- tree@ndlst
  tids <- tree@tips
  nids <- tree@nds
  rid <- tree@root
  # add first tip
  id <- tids[1]
  trstr <-  ''
  deja_vues <- rep(NA, length(ndlst))
  deja_vues[1] <- id
  spn <- ndlst[[id]][['spn']]
  dpth <- length(ndlst[[id]][['prid']])
  prid <- ndlst[[id]][['prid']]
  tpstr <- paste0(id, ':', spn)
  trstr <- paste0(rep('(', dpth), collapse='')
  trstr <- paste0(trstr, tpstr)
  # loop through nodes
  plyr::m_ply(2:(length(ndlst) - 1), .fun=tipBytip)
  trstr <- paste0(trstr, ');')
  write.table(x=trstr, file=file, quote=FALSE, row.names=FALSE,
              col.names=FALSE)
}

#' @name readTree
#' @title Read a Newick tree
#' @description Return a \code{TreeMan} or \code{TreeMen} object from a Newick treefile
#' @details Read a single or multiple trees from a file, or a text string. Parallelizable.
#' @param file file path
#' @param text Newick character string
#' @param ... \code{plyr} arguments
#' @seealso
#' \code{\link{writeTree}}, \code{\link{randTree}}, \url{https://en.wikipedia.org/wiki/Newick_format}
#' @useDynLib treeman
#' @useDynLib treeman findPrids
#' @export
#' @examples
#' library(treeman)
#' tree <- readTree(text="((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);")
readTree <- function(file=NULL, text=NULL, ...) {
  if(!is.null(file)) {
    trstr <- scan(file, what="raw", quiet=TRUE)
  } else {
    trstr <- text
  }
  if(length(trstr) > 1) {
    trees <- plyr::mlply(trstr, .fun=.readTree, ...)
    tree <- as(trees, 'TreeMen')
  } else {
    tree <- .readTree(trstr)
  }
  tree
}

.readTree <- function(trstr) {
  # Internals
  .idspn <- function(i) {
    mtdt <- substr(trstr, start=nds[i-1] + 1, stop=nds[i])
    mtdt <- gsub("(\\(|\\)|,)", "", mtdt)
    mtdt <- strsplit(mtdt, ":")[[1]]
    if(length(mtdt) == 0) {
      id <- paste0('n', i)
      spn <- NA
    } else if(mtdt[1] == ';') {
      id <- paste0('n', i)
      spn <- NA
    } else if(length(mtdt) == 1) {
      id <- mtdt
      spn <- NA
    } else if(length(mtdt) > 1 && mtdt[1] != "") {
      id <- mtdt[1]
      spn <- as.numeric(mtdt[2])
    } else {
      id <- paste0('n', i)
      spn <- as.numeric(mtdt[2])
    }
    c(id, spn)
  }
  # get nodes from string
  nds <- c(1, as.integer(gregexpr("(,|\\))", trstr)[[1]]) - 1)
  if(grepl(");", trstr)) {
    nds <- c(nds, nchar(trstr))
  }
  # get id and spn
  mtdt <- sapply(2:length(nds), FUN=.idspn)
  ids <- mtdt[1, ]
  spns <- as.numeric(mtdt[2, ])
  rm(mtdt)
  nds <- nds[-1]
  # gen prids
  opns <- gregexpr("\\(", trstr)[[1]]
  clss <- gregexpr("\\)", trstr)[[1]]
  prids <- .Call("findPrids", PACKAGE="treeman",
                 as.integer(nds),
                 as.integer(clss),
                 as.integer(opns))
  if(sum(prids == -1) > 1) {
    stop('Invalid tree string')
  }
  root <- which(prids == -1)
  prids <- match(prids, nds)
  tids <- which(!1:length(ids) %in% prids)
  prids[is.na(prids)] <- -1
  spns[is.na(spns)] <- -1
  treels <- list('ids'=ids, 'prids'=prids, 'spns'=spns,
                 'tids'=tids)
  rm(prids, tids, spns, ids)
  ndlst <- .updateNdlst(treels)
  # correct root
  if(length(root) > 0) {
    rid <- treels[['ids']][root]
    if(length(ndlst[[rid]][['spn']]) > 0) {
      ndlst[[rid]][['spn']] <- 0
      ndlst[[rid]][['pd']] <- sum(treels[['spns']])
    }
    ndlst[[rid]][['prid']] <- NULL
    tree <- new('TreeMan', ndlst=ndlst, root=rid)
  } else {
    tree <- new('TreeMan', ndlst=ndlst)
  }
  .updateTreeSlts(tree)
}
