
.getNdsPrdstsFrmMtrx <- function(ndmtrx, all_ids, ids, spns) {
  res <- apply(ndmtrx[ ,all_ids %in% ids],
               2, function(x) sum(spns[x]))
  res <- res + spns[all_ids %in% ids]
  names(res) <- ids
  res
}

.getNdsKidsFrmMtrx <- function(ndmtrx, all_ids, ids, tids) {
  res <- apply(ndmtrx[all_ids %in% ids, ],
               1, function(x) all_ids[x & all_ids %in% tids])
  names(res) <- ids
  res
}
