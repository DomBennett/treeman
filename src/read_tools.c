#include <R.h>
#include <Rinternals.h>

// Get kids
SEXP kids(SEXP nids_, SEXP tids_, SEXP prids_)
{
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
}

// Get prdst
SEXP prdst(SEXP nids_, SEXP prids_, SEXP spns_)
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
}

// Get pds
SEXP pd(SEXP nids_, SEXP prids_, SEXP spns_)
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
  int k;
  int j;
  int cc;
  int nxt_cc;
  int ids[nids];
  int nxt_ids[nids];
  for(i=0;i<(nids-1); i++) {
    double pd = 0;
    ids[0] = i;
    cc = 1;
    while(cc>0) {
      nxt_cc = 0;
      for(j=0;j<cc;j++) {
        for(k=0;k<(nids-1); k++) {
          if((prids[k]-1) == ids[j]) {
            pd += spns[j];
            nxt_ids[nxt_cc] = k;
            nxt_cc++;
          }
        }
      }
      cc = nxt_cc;
      for(j=0;j<cc;j++) {
        ids[j] = nxt_ids[j];
      }
    }
    REAL(res)[i] += pd;
  }
  UNPROTECT(1);
  return res;
}

// Get prids
SEXP prids(SEXP nds_, SEXP clss_, SEXP opns_)
{
  SEXP res;
  int n = length(nds_);
  int nclss = length(clss_);
  int nopns = length(opns_);
  int* nds = INTEGER(nds_);
  int* opns = INTEGER(opns_);
  int* clss = INTEGER(clss_);
  PROTECT(res=allocVector(REALSXP, n));
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
      REAL(res)[i] = -1;
    } else {
      for(j=0;j<n; j++) {
        if(nds[j] >= cls) {
          REAL(res)[i] = nds[j];
          break;
        }
      }
    }
  }
  UNPROTECT(1);
  return res;
}