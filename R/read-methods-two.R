library(inline)

# TODO:
# test with unrooted trees
# add more than just spns, nodelabels as well
# fix prid across treeman

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

.kids <- cfunction(c(nids_='integer', tids_='integer', prids_='integer'),
'
  SEXP res;
  int nids = asInteger(nids_);
  int ntips = length(tids_);
  int* tids = INTEGER(tids_);
  int* prids = INTEGER(prids_);
  PROTECT(res=allocMatrix(INTSXP, nids, ntips));
  int n = length(res);
  int i;
  for(i=0;i<n; i++) {
    INTEGER(res)[i] = 0;
  }
  for(i=0;i<ntips; i++) {
    int tid = tids[i] - 1;
    int id = prids[tid] - 1;
    while(id != nids-1) {
      INTEGER(res)[id + i * nids] = 1;
      id = prids[id] - 1;
    }
  }
  UNPROTECT(1);
  return res;
')

.prdst <- cfunction(c(nids_='integer', prids_='integer', spns_='numerical'),
'
  SEXP res;
  int nids = asInteger(nids_);
  double* spns = REAL(spns_);
  int* prids = INTEGER(prids_);
  PROTECT(res=allocVector(REALSXP, nids));
  int n = length(res);
  int i;
  for(i=0;i<n; i++) {
    REAL(res)[i] = 0;
  }
  for(i=0;i<(nids-1); i++) {
    double spn = spns[i];
    int id = prids[i] - 1;
    while(id != nids-1) {
      spn += spns[id];
      id = prids[id] - 1;
    }
    REAL(res)[i] += spn;
  }
  UNPROTECT(1);
  return res;
')

.pd <- function(ids, prids, spns) {
  res <- vector(length=length(ids))
  for(i in 1:(length(ids) - 1)) {
    ptids <- i
    pd <- 0
    while(length(ptids) > 0) {
      ptids <- which(prids %in% ptids)
      pd <- pd + sum(spns[ptids])
    }
    res[i] <- pd
  }
  res[length(ids)] <- sum(spns, na.rm=TRUE)
  res
}

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
  tids <- which(!1:length(ids) %in% prids)
  root <- length(prids)
  prids <- prids
  # generate other data
  kids <- .kids(as.integer(length(ids)), as.integer(tids), as.integer(prids[-root]))
  kids[root, ] <- 1
  kids <- kids == 1
  if(sum(is.na(spns)) == 1) {
    prdsts <- .prdst(as.integer(length(ids)), as.integer(prids[-root]),
                     as.numeric(spns[-root]))
    pds <- .pd(ids, prids, spns)
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
}