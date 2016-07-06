#include <R.h>
#include <Rinternals.h>

// Get prids
SEXP getPrids(SEXP prid_, SEXP prids_)
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