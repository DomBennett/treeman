.genMsg <- function(res, ids) {
  bad <- which(!res)
  msg <- ''
  for(i in bad[-length(res)]) {
    msg <- paste0(msg, ids[i], ', ')
  }
  paste0(msg, ids[bad[length(bad)]])
}

checkTree <- function(tree) {
  # TODO
  ids <- names(tree@ndlst)
  reported_ids <- as.vector(sapply(tree@ndlst, function(x) x[['id']]))
  cat("Checking IDs .... ")
  res <- ids == reported_ids
  if(!all(res)) {
    msg <- .genMsg(res, ids)
    cat('\n', msg, '\n')
  } else {
    cat('done\n')
  }
}