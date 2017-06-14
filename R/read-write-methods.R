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
#' \url{https://en.wikipedia.org/wiki/Newick_format},
#' \code{\link{readTree}}, \code{\link{randTree}},
#' \code{\link{readTrmn}}, \code{\link{writeTrmn}},
#' \code{\link{saveTreeMan}}, \code{\link{loadTreeMan}}
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
  dpth <- length(getNdPrids(tree, id))
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
#' @param wndmtrx T/F add node matrix? Default FALSE.
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \url{https://en.wikipedia.org/wiki/Newick_format},
#' \code{\link{addNdmtrx}}, \code{\link{writeTree}},
#' \code{\link{randTree}}, \code{\link{readTrmn}}, \code{\link{writeTrmn}},
#' \code{\link{saveTreeMan}}, \code{\link{loadTreeMan}}
#' @export
#' @examples
#' library(treeman)
#' tree <- readTree(text="((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);")
readTree <- function(file=NULL, text=NULL, wndmtrx=FALSE, parallel=FALSE,
                     progress='none') {
  if(!is.null(file)) {
    trstr <- scan(file, what="raw", quiet=TRUE)
  } else {
    trstr <- text
  }
  if(length(trstr) > 1) {
    trstr <- as.list(trstr)
    trees <- plyr::mlply(trstr, .fun=.readTree, wndmtrx=wndmtrx,
                         .progress=progress, .parallel=parallel)
    names(trees) <- NULL
    tree <- as(trees, 'TreeMen')
  } else {
    tree <- .readTree(trstr, wndmtrx)
  }
  tree
}

#' @useDynLib treeman
#' @useDynLib treeman cFindPrids
.readTree <- function(trstr, wndmtrx) {
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
              ndmtrx=NULL, wtxnyms=FALSE,
              prinds=prinds, tinds=tinds)
  tree <- updateSlts(tree)
  if(wndmtrx) {
    tree <- addNdmtrx(tree)
  }
  tree
}

#' @name writeTrmn
#' @title Write a .trmn tree
#' @description Write to disk a \code{TreeMan} or \code{TreeMan} object using the .trmn treefile
#' @details Write a tree(s) to file using the .trmn format.
#' It is faster to read and write tree files using treeman with the .trmn file format.
#' In addition it is possible to encode more information than possible with the
#' Newick, e.g. any taxonomic information and additional slot names added to 
#' the tree are recorded in the file.
#' @param tree TreeMan object or TreeMen object
#' @param file file path
#' @seealso
#' \code{\link{readTrmn}},
#' \code{\link{readTree}},\code{\link{writeTree}},
#' \code{\link{randTree}}, \code{\link{saveTreeMan}}, \code{\link{loadTreeMan}}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' writeTrmn(tree, file='test.trmn')
#' tree <- readTrmn('test.trmn')
#' file.remove('test.trmn')
writeTrmn <- function(tree, file) {
  .unpack <- function(ntree) {
    .makeDataFrame(ntree, tree@treelst[[ntree]])
  }
  .makeDataFrame <- function(ntree, tree) {
    res <- data.frame(tree=ntree, prind=tree@prinds)
    res[['id']] <- names(tree@ndlst)
    if(tree@wspn) {
      res[['spn']] <- sapply(tree@ndlst, function(x) x[['spn']])
    }
    if(tree@wtxnyms) {
      res[['txnym']] <- sapply(tree@ndlst,
                               function(x) paste0(x[['txnym']], collapse='|'))
    }
    # add any additional slots
    if(length(tree@othr_slt_nms) > 0) {
      for(slt_nm in tree@othr_slt_nms) {
        res[[slt_nm]] <- sapply(tree@ndlst,
                                function(x) x[[slt_nm]])
      }
    }
    res
  }
  if('TreeMan' %in% is(tree)) {
    res <- .makeDataFrame(1, tree)
  } else if('TreeMen' %in% is(tree)) {
    res <- plyr::mdply(.data=data.frame(ntree=1:tree@ntrees),
                       .fun=.unpack)
    res <- res[ ,-1]
  } else {
    stop("`tree` must be TreeMan or TreeMen object.")
  }
  write.csv(res, file=file, quote=FALSE, row.names=FALSE)
}

#' @name readTrmn
#' @title Read a .trmn tree
#' @description Return a \code{TreeMan} or \code{TreeMen} object from a .trmn treefile
#' @details Read a tree(s) from a file using the .trmn format.
#' It is faster to read and write tree files using treeman with the .trmn file format.
#' In addition it is possible to encode more information than possible with the
#' Newick, e.g. any taxonomic information and additional slot names added to 
#' the tree are recorded in the file.
#' @param file file path
#' @param wndmtrx T/F add node matrix? Default FALSE.
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{writeTrmn}},
#' \code{\link{readTree}},\code{\link{writeTree}},
#' \code{\link{randTree}}, \code{\link{saveTreeMan}}, \code{\link{loadTreeMan}}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' writeTrmn(tree, file='test.trmn')
#' tree <- readTrmn('test.trmn')
#' file.remove('test.trmn')
readTrmn <- function(file, wndmtrx=FALSE, parallel=FALSE,
                     progress='none') {
  .pack <- function(i) {
    .readTrmn(inpt[inpt[['tree']] == i, ],
              wndmtrx)
    
  }
  inpt <- read.csv(file, stringsAsFactors=FALSE)
  trids <- unique(inpt[['tree']])
  trees <- plyr::mlply(.data=trids, .fun=.pack,
                       .parallel=parallel, .progress=progress)
  if(length(trees) == 1) {
    res <- trees[[1]]
  } else {
    trees <- trees[1:length(trees)]
    names(trees) <- NULL
    res <- as(trees, 'TreeMen')
  }
  res
}

.readTrmn <- function(inpt, wndmtrx) {
  .add <- function(i) {
    nd <- vector("list", length=4)
    names(nd) <- c('id', 'ptid', 'prid', 'spn')
    nd[['id']] <- ids[i]
    nd[['spn']] <- spns[i]
    nd[['prid']] <- ids[prinds[i]]
    nd[['ptid']] <- ptids[ptnds_pool == i]
    nd
  }
  prinds <- inpt[['prind']]
  # all internal nodes should occur more than once (twice for bifurcating trees)
  prind_test <- sum(prinds == 1:length(prinds)) == 1
  prind_test <- all(table(prinds) > 1) & prind_test
  if(!prind_test) {
    stop('Tree is corrupted, check node structure is hierarchical.')
  }
  ids <- inpt[['id']]
  if('spn' %in% names(inpt) && !is.na(inpt[['spn']][1])) {
    spns <- inpt[['spn']]
  } else {
    spns <- rep(0 , length(ids))
  }
  tinds <- which(!1:length(ids) %in% prinds)
  root <- which(1:length(prinds) == prinds)
  ptids <- ids[-root]
  ptnds_pool <- prinds[-root]
  ndlst <- lapply(1:length(ids), .add)
  names(ndlst) <- ids
  tree <- new('TreeMan', ndlst=ndlst, root=ids[root],
              ndmtrx=NULL, wtxnyms=FALSE,
              prinds=prinds, tinds=tinds)
  if('txnym' %in% names(inpt) && !is.na(inpt[['txnym']][1])) {
    txnyms <- strsplit(inpt[['txnym']], '\\|')
    names(txnyms) <- ids
    tree <- setTxnyms(tree, txnyms)
  }
  othr_slt_nms <- names(inpt)[!names(inpt) %in%
                                c('id', 'prind', 'spn', 'txnym', 'tree')]
  if(length(othr_slt_nms) > 0) {
    for(slt_nm in othr_slt_nms) {
      tree <- setNdsOther(tree, ids=inpt[['id']],
                          vals=inpt[[slt_nm]], slt_nm=slt_nm)
    }
  }
  tree <- updateSlts(tree)
  if(wndmtrx) {
    tree <- addNdmtrx(tree)
  }
  tree
}

#' @name saveTreeMan
#' @title Save a TreeMan object in serialization format
#' @description \code{TreeMan} equivalent to \code{save()} but able to handle
#' node matrices.
#' @details It is not possible to use \code{save()} on \code{TreeMan} objects
#' with node matrices. Node matrices are bigmemory matrices and are therefore outside
#' the R environment, see bigmemory documentation for more information. Saving and loading
#' a bigmemory matrix may cause memory issues in R and cause R to crash.
#' 
#' This function can safely store a \code{TreeMan} object with and without
#' a node matrix. This function stores the tree using the serialization format and the node
#' matrix as a hidden .csv. Both parts of the tree can be reloaded to an R environment
#' with \code{loadTreeMan()}. The hidden node matrix filename is based on the file argument:
#' \code{file + _ndmtrx}
#' 
#' Reading and writing trees with \code{saveTreeMan()} and
#' \code{loadTreeMan} is faster than any of the other read and write functions.
#' @param tree \code{TreeMan} object
#' @param file file path
#' @seealso
#' \code{\link{loadTreeMan}},
#' \code{\link{readTree}},\code{\link{writeTree}},
#' \code{\link{readTrmn}}, \code{\link{writeTrmn}}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(100, wndmtrx=TRUE)
#' saveTreeMan(tree, file='test.RData')
#' rm(tree)
#' tree <- loadTreeMan(file='test.RData')
#' file.remove('test.RData', 'testRData_ndmtrx')
saveTreeMan <- function(tree, file) {
  ndmtrx_file <- paste0(gsub('\\.', '', file), '_ndmtrx')
  if(!is.null(tree@ndmtrx)) {
    bigmemory::write.big.matrix(x=tree@ndmtrx, filename=ndmtrx_file)
    tree <- rmNdmtrx(tree)
  }
  save(list=c('tree', 'ndmtrx_file'), file=file)
}

#' @name loadTreeMan
#' @title Load a TreeMan object in serialization format
#' @description \code{TreeMan} equivalent to \code{load()} but able to handle
#' node matrices.
#' @details It is not possible to use \code{save()} on \code{TreeMan} objects
#' with node matrices. Node matrices are bigmemory matrices and are therefore outside
#' the R environment, see bigmemory documentation for more information. Saving and loading
#' a bigmemory matrix may cause memory issues in R and cause R to crash.
#' 
#' This function can safely read a \code{TreeMan} object with and without
#' a node matrix. \code{saveTreeMan()} function stores the tree using the serialization format
#' and the node matrix as a hidden .csv. Both parts of the tree can be reloaded to an R environment
#' with \code{loadTreeMan()}. The hidden node matrix filename is based on the file argument:
#' \code{file + _ndmtrx}
#' 
#' Reading and writing trees with \code{saveTreeMan()} and
#' \code{loadTreeMan} is faster than any of the other read and write functions.
#' @param file file path
#' @seealso
#' \code{\link{saveTreeMan}},
#' \code{\link{readTree}},\code{\link{writeTree}},
#' \code{\link{readTrmn}}, \code{\link{writeTrmn}}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(100, wndmtrx=TRUE)
#' saveTreeMan(tree, file='test.RData')
#' rm(tree)
#' tree <- loadTreeMan(file='test.RData')
#' file.remove('test.RData', 'testRData_ndmtrx')
loadTreeMan <- function(file) {
  ndmtrx_file <- NULL
  load(file)
  if(file.exists(ndmtrx_file)) {
    tree@ndmtrx <- bigmemory::read.big.matrix(filename=ndmtrx_file,
                                              type='integer', shared=FALSE)
  }
  tree
}