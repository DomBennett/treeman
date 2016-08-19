#include <R.h>
#include <Rinternals.h>

// Return vector of position of prid in Newick string
SEXP cFindPrids(SEXP nds_, SEXP clss_, SEXP opns_)
{
  SEXP res_prids;
  int n = length(nds_);
  int nclss = length(clss_);
  int nopns = length(opns_);
  int* nds = INTEGER(nds_);
  int* opns = INTEGER(opns_);
  int* clss = INTEGER(clss_);
  PROTECT(res_prids=allocVector(REALSXP, n));
  int i, j, cls_cc, opn_cc, cls;
  int cls_up[nclss];
  int opn_up[nopns];
  for(i=0;i<n; i++) {
    //find preceding node by finding closing bracket
    cls_cc = 0;
    opn_cc = 0;
    cls = -1;
    //find all closing brackets upstream of node
    for(j=0;j<nclss; j++) {
      if(nds[i] < clss[j]) {
        cls_up[cls_cc] = clss[j];
        cls_cc++;
      }
    }
    //find all opening brackets upstream of node
    for(j=0;j<nopns; j++) {
      if(nds[i] < opns[j]) {
        opn_up[opn_cc] = opns[j];
        opn_cc++;
      }
    }
    for(j=0;j<cls_cc; j++) {
      if(opn_cc < j + 1) {
        cls = cls_up[j];
        break;
      }
      if(opn_up[j] > cls_up[j]) {
        cls = cls_up[j];
        break;
      }
    }
    if(cls == -1) {
      //use -1 for NA
      REAL(res_prids)[i] = -1;
    } else {
      for(j=0;j<n; j++) {
        if(nds[j] >= cls) {
          REAL(res_prids)[i] = nds[j];
          break;
        }
      }
    }
  }
  UNPROTECT(1);
  return res_prids;
}