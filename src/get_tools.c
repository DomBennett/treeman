#include <R.h>
#include <Rinternals.h>

// Loops from qry ids to root
// Return matrix of 01s for presence/absence
// all nodes are rows, qry ids are cols
SEXP cGetNdmtrx(SEXP nids_, SEXP qrys_, SEXP prids_)
{
  SEXP ndmtrx;
  int nids = asInteger(nids_);
  int nqrys = length(qrys_);
  int* qrys = INTEGER(qrys_);
  int* prids = INTEGER(prids_);
  PROTECT(ndmtrx=allocMatrix(INTSXP, nids, nqrys));
  int n = length(ndmtrx);
  int i;
  for(i=0;i<n; i++) {
    INTEGER(ndmtrx)[i] = 0;
  }
  int qry;
  int id;
  int prv_ids[2];
  for(i=0;i<nqrys; i++) {
    qry = qrys[i] - 1;
    id = prids[qry] - 1;
    // stop while loop by testing for self-reference
    // if the last two assignments to id are the same, break
    prv_ids[0] = -1;
    prv_ids[1] = id;
    while(prv_ids[0] != prv_ids[1]) {
      prv_ids[0] = id;
      INTEGER(ndmtrx)[id + i * nids] = 1;
      id = prids[id] - 1;
      prv_ids[1] = id;
    }
  }
  UNPROTECT(1);
  return ndmtrx;
}

// Get prids for a node
SEXP cGetNdPrids(SEXP prid_, SEXP prids_)
{
  int nprids = length(prids_);
  int init_res[nprids+2];
  int prid = asInteger(prid_);
  //vector of internal node prids
  // duplicated and protected
  int* prids = INTEGER(PROTECT(duplicate(prids_)));
  int i=2;
  init_res[0] = -1;
  init_res[1] = -1;
  //loop through until either of
  //the previous two prids equal current prid
  while(init_res[i-1] != prid &
        init_res[i-2] != prid) {
    init_res[i] = prid;
    prid = prids[prid - 1];
    i++;
  }
  SEXP res;
  int j;
  PROTECT(res=allocVector(INTSXP, i-2));
  for(j=0;j<(i-2);j++) {
    INTEGER(res)[j] = init_res[j+2];
  }
  UNPROTECT(1);
  return res;
}

// get ptids for a node
SEXP cGetNdPtids(SEXP id_, SEXP prids_) {
  SEXP res_ptids;
  R_xlen_t nids = xlength(prids_);
  int id = asInteger(id_);
  // vector of internal node prids
  // duplicate, potential to be a C generated object via cFindPrids
  // without duplicating, this function will modify other vectors
  // in the R environment
  // Protect duplicate again, no longer argument
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
}