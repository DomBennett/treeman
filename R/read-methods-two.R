# Test strings
trstr <- scan(file="/Users/djb208/Coding/Project-EPI/0_data/trees/mammals.tre", what='raw')
trstr <- "((A:1.0,(A1:1.0,A2:1.0):1.0):1.0,(C:1.0,D:1.0):1.0);"

# Internal functions for conversion to C
.prids <- function(trstr) {
  # get prids
  prids <- rep(NA, length(nds))
  opns <- gregexpr("\\(", trstr)[[1]]
  clss <- gregexpr("\\)", trstr)[[1]]
  for(i in 1:length(nds)) {
    # convert to C
    # find preceding node by finding closing bracket
    cls_dwn <- clss[nds[i] < clss]
    opn_dwn <- opns[nds[i] < opns]
    for(j in 1:length(cls_dwn)) {
      if(length(opn_dwn) < j) {
        break
      }
      if(opn_dwn[j] > cls_dwn[j]) {
        break
      }
    }
    prids[i] <- nds[cls_dwn[j] < nds][1]
  }
  prids
}

.mtdt <- function(nds) {
  res <- matrix(FALSE, nrow=length(nds), ncol=2)
  for(i in 1:length(nds)) {
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
    res[i, 1] <- id
    res[i, 2] <- spn
  }
  res
}

.kids <- function(tids, prids) {
  nids <- unique(prids[!is.na(prids)])
  # loop through tids, step through nids with while loop
  # whenever a nid is encountered, select TRUE for that nid
  res <- matrix(FALSE, nrow=length(nids),
                ncol=length(tids))
  for(i in 1:length(tids)) {
    nid <- prids[ids == tids[i]]
    while(TRUE) {
      res[nid == nids, i] <- TRUE
      nid <- prids[ids == nid]
      if(is.na(nid)) {
        break
      }
    }
  }
  res
}

.prdst <- function(ids, prids, spns, root) {
  res <- rep(0, length(ids))
  for(i in 1:(length(prids) - 1)) {
    res[i] <- res[i] + spns[ids == ids[i]]
    id <- prids[i]
    while(TRUE) {
      if(id == root) {
        break
      }
      res[i] <- res[i] + spns[ids == id]
      id <- prids[ids == id]
    }
  }
  res
}

.pd <- function() {
  cat("TODO")
}

.addwospn <- function(i) {
  nd <- vector("list", length=4)
  nd[['id']] <- ids[i]
  nd[['kids']] <- tids[kids[nids == ids[i], ]]
  nd[['prid']] <- prids[i]
  ptids <- ids[prids == ids[i]]
  nd[['ptid']] <- ptids[!is.na(ptids)]
  nd
}

.addwspn <- function(i) {
  nd <- vector("list", length=7)
  nd[['id']] <- ids[i]
  nd[['spn']] <- spns[i]
  nd[['prdst']] <- prdsts[i]
  nd[['pd']] <- pds[i]
  nd[['kids']] <- tids[kids[nids == ids[i], ]]
  nd[['prid']] <- prids[i]
  ptids <- ids[prids == ids[i]]
  nd[['ptid']] <- ptids[!is.na(ptids)]
  nd
}

.read <- function(trstr) {
  nds <- gregexpr(":", trstr)[[1]]
  if(grepl(");", trstr)) {
    nds <- c(nds, nchar(trstr))
  }
  # get id and spn
  nonmtdt <- gregexpr("(\\(|\\)|,)", trstr)[[1]]
  mtdt <- .mtdt(nds)
  ids <- mtdt[ ,1]
  spns <- as.numeric(mtdt[ ,2])
  # gen prids
  prids <- .prids(trstr)
  prids <- match(prids, nds)
  prids <- ids[prids]
  tids <- ids[!ids %in% prids]
  root <- ids[length(prids)]
  # generate other data
  kids <- .kids(tids, prids)
  if(is.na(spns) == 1) {
    prdsts <- .prdst(ids, prids, spns, root)
    pds <- .pd()
  }
  # generate ndlst
  if(is.na(spns) == 1) {
    ndlst <- lapply(1:length(ids), .addwospn)
  } else {
    ndlst <- lapply(1:length(ids), .addwspn)
  }
  # create TreeMan object
  tree <- new('TreeMan', ndlst=ndlst, root=root)
  .updateTreeSlts(tree)
}