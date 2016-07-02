library(inline)

# TODO:
# test with unrooted trees
# add more than just spns, nodelabels as well
# fix prid across treeman

# Test strings
trstr <- "((A:1.0,(A1:1.0,A2:1.0):1.0):1.0,(C:1.0,D:1.0):1.0);"
trstr <- scan(file="/Users/djb208/Coding/Project-EPI/0_data/trees/mammals.tre", what='raw')

.read <- function(trstr) {
  # Internals
  .addwospn <- function(i) {
    nd <- vector("list", length=4)
    names(nd) <- c('id', 'ptid', 'prid', 'kids',
                   'spn', 'pd', 'prdst')
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
  .idspn <- function(i) {
    uppr <- which(nonmtdt < nds[i])
    uppr <- uppr[length(uppr)]
    lwr <- which(nonmtdt > nds[i])[1]
    mtdt <- substr(trstr, start=nonmtdt[uppr] + 1,
                   stop=nonmtdt[lwr] - 1)
    mtdt <- strsplit(mtdt, ":")[[1]]
    if(length(mtdt) > 1 && mtdt[1] != "") {
      id <- mtdt[1]
      spn <- as.numeric(mtdt[2])
    } else {
      id <- paste0('n', i)
      spn <- as.numeric(mtdt[2])
    }
    c(id, spn)
  }
  # get nodes from string
  nds <- gregexpr(":", trstr)[[1]]
  if(grepl(");", trstr)) {
    nds <- c(nds, nchar(trstr))
  }
  # get id and spn
  nonmtdt <- gregexpr("(\\(|\\)|,)", trstr)[[1]]
  mtdt <- sapply(1:length(nds), FUN=.idspn)
  ids <- mtdt[1, ]
  spns <- as.numeric(mtdt[2, ])
  # gen prids
  opns <- gregexpr("\\(", trstr)[[1]]
  clss <- gregexpr("\\)", trstr)[[1]]
  prids <- .Call("prids", PACKAGE="treeman",
                 as.integer(nds),
                 as.integer(clss),
                 as.integer(opns))
  if(sum(prids == -1) > 1) {
    stop('Invalid tree string')
  }
  prids <- match(prids, nds)
  tids <- which(!1:length(ids) %in% prids)
  root <- length(prids)
  prids <- prids
  # generate other data
  kids <- .Call("kids", PACKAGE="treeman",
                as.integer(length(ids)),
                as.integer(tids),
                as.integer(prids[-root]))
  kids[root, ] <- 1
  kids <- kids == 1
  if(sum(is.na(spns)) == 1) {
    prdsts <- .Call("prdst", PACKAGE="treeman",
                    as.integer(length(ids)),
                    as.integer(prids[-root]),
                    as.numeric(spns[-root]))
    pds <- .Call("pd", PACKAGE="treeman",
                 as.integer(length(ids)),
                 as.integer(prids[-root]),
                 as.numeric(spns[-root]))
    pds[root] <- sum(spns[-root])
  }
  # generate ndlst
  tids <- ids[tids]
  prids <- ids[prids]
  if(sum(is.na(spns)) > 1) {
    ndlst <- lapply(1:length(ids), .addwospn)
  } else {
    ndlst <- lapply(1:length(ids), .addwspn)
  }
  names(ndlst) <- ids
  root <- ids[root]
  ndlst[[root]][['spn']] <- 0
  ndlst[[root]][['prid']] <- NULL
  # create TreeMan object
  tree <- new('TreeMan', ndlst=ndlst, root=root)
  tree <- treeman:::.updateTreeSlts(tree)
  tree
}

readTree <- function(file=NULL, text=NULL, ...) {
  if(!is.null(file)) {
    trstr <- scan(file, what="raw", quiet=TRUE)
  } else {
    trstr <- text
  }
  if(length(trstr) > 1) {
    trees <- plyr::mlply(trstr, .fun=.read, ...)
    tree <- as(trees, 'TreeMen')
  } else {
    tree <- .read(trstr)
  }
  tree
}