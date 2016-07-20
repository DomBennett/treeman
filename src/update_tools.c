#include <R.h>
#include <Rinternals.h>

// Return vector of prdsts of every node
// Only returns true numbers if nds vector provided are rootable
SEXP getPrdstVec(SEXP nids_, SEXP prids_, SEXP spns_)
{
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
  for(i=0;i<nids; i++) {
    double spn = spns[i];
    int id = prids[i] - 1;
    while(id != -2) {
      spn += spns[id];
      id = prids[id] - 1;
    }
    REAL(res)[i] += spn;
  }
  UNPROTECT(1);
  return res;
}

// Return vector of pds for every node
// Only returns true numbers if nds vector provided are rootable
SEXP getPdVec(SEXP nids_, SEXP prids_, SEXP spns_)
{
  int nids = asInteger(nids_);
  double* spns = REAL(spns_);
  int* prids = INTEGER(prids_);
  int i;
  double res[nids];
  for(i=0;i<nids; i++) {
    res[i] = 0;
  }
  for(i=0;i<nids; i++) {
    double spn = spns[i];
    int id = prids[i] - 1;
    while(id != -2) {
      res[id] += spn;
      id = prids[id] - 1;
    }
  }
  SEXP out;
  PROTECT(out=allocVector(REALSXP, nids));
  for(i=0;i<nids;i++) {
    REAL(out)[i] += res[i];
  }
  UNPROTECT(1);
  return out;
}