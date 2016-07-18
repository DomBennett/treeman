
updateTree <- function(tree) {
  # Update the slots for a tree
  wo_pstndes <- sapply(tree@ndlst,
                       function(n) length(n[['ptid']]) == 0)
  tree@tips <- sort(names(wo_pstndes)[wo_pstndes])
  tree@ntips <- length(tree@tips)
  tree@nds <- sort(names(wo_pstndes)[!wo_pstndes])
  tree@nnds <- length(tree@nds)
  tree@all <- c(tree@tips, tree@nds)
  tree@nall <- length(tree@all)
  if(length(tree@root) > 0) {
    wspn <- names(tree@ndlst)[names(tree@ndlst) != tree@root]
  } else {
    wspn <- names(tree@ndlst)
  }
  tree@wspn <- all(sapply(tree@ndlst[wspn], function(n) !is.null(n[['spn']])))
  if(tree@wspn) {
    if(length(tree@root) > 0) {
      tip_prdsts <- sapply(tree@tips, .getPrdst, ndlst=tree@ndlst)
      tree@age <- max(tip_prdsts)
      extant_is <- (tree@age - tip_prdsts) <= tree@tol
      tree@ext <- names(extant_is)[extant_is]
      tree@exc <- tree@tips[!tree@tips %in% tree@ext]
      tree@ultr <- all(tree@tips %in% tree@ext)
    } else {
      tree@ext <- tree@exc <- vector()
      tree@ultr <- FALSE
      tree@age <- numeric()
    }
    tree@pd <- sum(sapply(tree@ndlst, function(n) n[['spn']]))
  } else {
    tree@age <- tree@pd <- numeric()
    tree@ext <- tree@ext <- vector()
    tree@ultr <- logical()
  }
  tree@ply <- any(sapply(tree@ndlst, function(n) length(n[['ptid']]) > 2))
  initialize(tree)
}