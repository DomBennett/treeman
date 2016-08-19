#' @name writeTree
#' @title Write a Newick tree
#' @description Creates a Newick tree from a \code{TreeMan} object.
#' @details The \code{ndLabels} argument can be used to add a user defined node label in
#' the Newick tree. It should take only 1 argument, \code{nd}, the node represented as a list.
#' It should only return a single character value that can be added to a newick string.
#' @param tree \code{TreeMan} object
#' @param file file path
#' @param append T/F append tree to already existing file
#' @param ndLabels node label function
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{readTree}}, \code{\link{randTree}}, \url{https://en.wikipedia.org/wiki/Newick_format}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' ndLabels <- function(n) {
#' paste0(n[['id']], '_ndlabel')
#' }
#' writeTree(tree, file='example.tre', ndLabels=ndLabels)
#' file.remove('example.tre')
writeTree <- function(tree, file, append=FALSE, ndLabels=function(nd){
  return(NULL)
  }, parallel=FALSE, progress="none") {
  if(is(tree) == 'TreeMen') {
    plyr::m_ply(tree@treelst, .fun=.writeTree, file=file,
                append=TRUE, ndLabels=ndLabels,
                .progress=progress, .parallel=parallel)
  } else if(is(tree) == "TreeMan") {
    .writeTree(tree, file, append, ndLabels)
  } else {
    stop('`tree` must be TreeMan or TreeMen')
  }
  NULL
}

.writeTree <- function(tree, file, append, ndLabels) {
  tipBytip <- function(i) {
    kids <- getNdKids(tree, prid)
    ids <- c(kids, prid, ndlst[[prid]][['prid']])
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
              col.names=FALSE, append=append)
}

#' @name readTree
#' @title Read a Newick tree
#' @description Return a \code{TreeMan} or \code{TreeMen} object from a Newick treefile
#' @details Read a single or multiple trees from a file, or a text string. Parallelizable
#' when reading multiple trees.
#' @param file file path
#' @param text Newick character string
#' @param update T/F update tree slots after generation? Default TRUE.
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{writeTree}}, \code{\link{randTree}}, \url{https://en.wikipedia.org/wiki/Newick_format}
#' @export
#' @examples
#' library(treeman)
#' tree <- readTree(text="((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);")
readTree <- function(file=NULL, text=NULL, update=TRUE, parallel=FALSE,
                     progress='none') {
  if(!is.null(file)) {
    trstr <- scan(file, what="raw", quiet=TRUE)
  } else {
    trstr <- text
  }
  if(length(trstr) > 1) {
    trstr <- as.list(trstr)
    trees <- plyr::mlply(trstr, .fun=.readTree, update=update,
                         .progress=progress, .parallel=parallel)
    tree <- as(trees, 'TreeMen')
  } else {
    tree <- .readTree(trstr, update)
  }
  tree
}

#' @useDynLib treeman
#' @useDynLib treeman cFindPrids
.readTree <- function(trstr, update) {
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
  .add <- function(i) {
    nd <- vector("list", length=4)
    names(nd) <- c('id', 'ptid', 'prid', 'spn')
    nd[['id']] <- ids[i]
    nd[['spn']] <- spns[i]
    nd[['prid']] <- ids[prinds[i]]
    nd[['ptid']] <- ptids[ptnds_pool == i]
    nd
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
  prinds <- .Call("cFindPrids", PACKAGE="treeman",
                  as.integer(nds),
                  as.integer(clss),
                  as.integer(opns))
  if(sum(prinds == -1) > 1) {
    stop('Invalid tree string')
  }
  root <- which(prinds == -1)
  prinds <- match(prinds, nds)
  tinds <- which(!1:length(ids) %in% prinds)
  prinds[is.na(prinds)] <- root
  spns[is.na(spns)] <- 0
  ptids <- ids[-root]
  ptnds_pool <- prinds[-root]
  ndlst <- lapply(1:length(ids), .add)
  names(ndlst) <- ids
  tree <- new('TreeMan', ndlst=ndlst, root=ids[root],
              ndmtrx=bigmemory::big.matrix(1,1),
              prinds=prinds, tinds=tinds)
  if(update) {
    tree <- updateTree(tree)
  } else {
    # init basic slots
    tree@updtd <- FALSE
    tree@tips <- sort(ids[tinds])
    tree@ntips <- length(tinds)
    tree@nds <- sort(ids[ids != tree@tips])
    tree@nnds <- length(tree@nds)
    tree@all <- names(tree@ndlst)
    tree@nall <- length(tree@all)
    tree@wspn <- any(spns > 0)
  }
  tree
}




# TODO: develop .trmn file format