# TODO:
# recreating seg fault error

library(treeman)
library(inline)


cGetNdPtids<- cfunction(c('id_'='integer', 'prids_'='integer'), {
  '
  SEXP res_ptids;
  R_xlen_t nids = xlength(prids_);
  int id = asInteger(id_);
  int* prids = INTEGER(PROTECT(duplicate(prids_)));
  PROTECT(res_ptids=allocVector(INTSXP, nids));
  int *pres = INTEGER(res_ptids);
  int qrys[nids+1];
  // init res_ptids and qrys
  memset(pres, 0, nids * sizeof(int));
  memset(qrys, -1, (nids + 1) * sizeof(int));
  int qry=id;
  int ni=0;
  int nqrys=0;
  R_xlen_t i;
  while(qry != -1) {
  // remove qry from prids
  for(i=0;i<nids; i++) {
  if(qry == (i + 1)) {
  prids[i] = -1;
  }
  }
  // search for qry in prids
  for(i=0;i<nids; i++) {
  if(qry == prids[i]) {
  pres[i] = 1;
  nqrys = nqrys + 1;
  qrys[nqrys] = i + 1;
  }
  }
  // update qry
  ni = ni + 1;
  qry = qrys[ni];
  }
  UNPROTECT(2);
  return res_ptids;
  '
})

getPtids <- function(ids, nids, prinds) {
  calc <- function(id) {
    res <- cGetNdPtids(as.integer(id), as.integer(prinds))
    nids[which(res > 0)]
  }
  l_data <- data.frame(id=ids, stringsAsFactors=FALSE)
  out <- plyr::mlply(.data=l_data, .fun=calc)
  names(out) <- attr(out, 'split_labels')[,1]
  res <- out[1:length(out)]
  res
}

options(CBoundsCheck = TRUE)
gctorture(FALSE)

tree <- randTree(1000, update=FALSE)
ndlst <- tree@ndlst
prinds <- tree@prinds
nids <- names(ndlst)
ids <- which(tree@all %in% tree@nds)

res_1 <- getPtids(ids, nids, prinds)
res_2 <- getPtids(ids, nids, prinds)

for(i in 1:length(res_1)) {
  if(!(all(res_1[[i]] %in% res_2[[i]]) &
       all(res_2[[i]] %in% res_1[[i]]))) {
    print(i)
  }
}


any(sapply(res_1, length) == 0)
any(sapply(res_2, length) == 0)
any(sapply(res_3, length) == 0)
which(sapply(res_1, length) == 0)
which(sapply(res_2, length) == 0)
which(sapply(res_3, length) == 0)

