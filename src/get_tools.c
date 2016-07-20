#include <R.h>
#include <Rinternals.h>

// Loops from qry ids to root
// Return matrix of 01s for presence/absence
// all nodes are rows, qry ids are cols
SEXP cGetNdsMat(SEXP nids_, SEXP qrys_, SEXP prids_)
{
  SEXP res;
  int nids = asInteger(nids_);
  int nqrys = length(qrys_);
  int* qrys = INTEGER(qrys_);
  int* prids = INTEGER(prids_);
  PROTECT(res=allocMatrix(INTSXP, nids, nqrys));
  int n = length(res);
  int i;
  for(i=0;i<n; i++) {
    INTEGER(res)[i] = 0;
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
      INTEGER(res)[id + i * nids] = 1;
      id = prids[id] - 1;
      prv_ids[1] = id;
    }
  }
  UNPROTECT(1);
  return res;
}

// Get prids for a node
SEXP cGetNdPrids(SEXP prid_, SEXP prids_)
{
  int nprids = length(prids_);
  int init_res[nprids+2];
  int prid = asInteger(prid_);
  int* prids = INTEGER(prids_);  //vector of internal node prids
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
  int nids = length(prids_);
  int id = asInteger(id_);
  int* prids = INTEGER(prids_);  //vector of internal node prids
  SEXP res;
  PROTECT(res=allocVector(INTSXP, nids));
  int qrys[nids+1];
  int i;
  // init res and qrys
  for(i=0;i<nids; i++) {
    INTEGER(res)[i] = 0;
    qrys[i] = -1;
  }
  qrys[nids+1] = -1;
  int qry=id;
  int ni=0;
  int nqrys=0;
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
        INTEGER(res)[i] = 1;
        nqrys = nqrys + 1;
        qrys[nqrys] = i + 1;
      }
    }
    // update qry
    ni = ni + 1;
    qry = qrys[ni];
  }
  UNPROTECT(1);
  return res;
}