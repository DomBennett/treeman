#TODO:
# .updatePD
# .updatePrdst
# .updateChildren
# .updateSpnData
# .updateAll

.updateSlots <- function(tree) {
  wo_pstndes <- sapply(tree@nodelist,
                       function(n) length(n[['ptid']]) == 0)
  tree@tips <- sort(names(wo_pstndes)[wo_pstndes])
  tree@ntips <- length(tree@tips)
  tree@nodes <- sort(names(wo_pstndes)[!wo_pstndes])
  tree@nnodes <- length(tree@nodes)
  tree@all <- c(tree@tips, tree@nodes)
  tree@nall <- length(tree@all)
  wspn <- names(tree@nodelist)[names(tree@nodelist) != tree@root]
  tree@wspn <- all(sapply(tree@nodelist[wspn], function(n) !is.null(n[['span']])))
  if(tree@wspn) {
    if(!is.null(tree@root)) {
      tree@age <- max(sapply(tree@nodelist[wspn], function(n) n[['prdst']]))
      extant_is <- unlist(sapply(tree@tips, function(i) {
        (tree@age - tree@nodelist[[i]][['prdst']]) <= tree@tol}))
      tree@ext <- names(extant_is)[extant_is]
      tree@exc <- tree@tips[!tree@tips %in% tree@ext]
      tree@ultr <- all(tree@tips %in% tree@ext)
    }
    tree@pd <- tree@nodelist[[tree@root]][['pd']]
  } else {
    tree@age <- tree@pd <- numeric()
    tree@ext <- tree@ext <- vector()
    tree@ultr <- logical()
  }
  tree@ply <- any(sapply(tree@nodelist, function(n) length(n[['ptid']]) > 2))
  initialize(tree)
}