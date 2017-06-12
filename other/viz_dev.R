calcRootYs <- function(y, ptids) {
  # determine the possible ys for the ptids of the root
  hlfwy <- round(length(ptids)/2)
  pstvs <- sort(cumsum(dpths[ptids][1:hlfwy]),
                decreasing=TRUE)
  ngtvs <- cumsum(dpths[ptids][(hlfwy+1):length(ptids)])
  pstvs <- as.list(pstvs + y)
  strt <- y + 1
  for(i in length(pstvs):1) {
    pstvs[[i]] <- strt:pstvs[[i]]
    strt <- max(pstvs[[i]]) + 1
  }
  ngtvs <- as.list(y - ngtvs)
  strt <- y - 1
  for(i in 1:length(ngtvs)) {
    ngtvs[[i]] <- strt:ngtvs[[i]]
    strt <- min(ngtvs[[i]]) - 1
  }
  c(pstvs, ngtvs)
}
calcYs <- function(y2, ys, id, ptids) {
  pssbls <- list(max(ys[[id]]):(y2+1), min(ys[[id]]):(y2-1))
  ns <- sapply(pssbls, length)
  for(ptid in ptids) {
    dpths[[ptid]]
  }
  
}
calcY <- function(ys, ptids) {
  res <- rep(NA, length(ptids))
  for(i in 1:length(ptids)) {
    ptid <- ptids[[i]]
    ptid_dpth <- dpths[[ptid]]
    if(ptid_dpth == 1) {
      res[i] <- ys[[ptid]]
    } else {
      pptid <- tree@ndlst[[ptid]][['ptid']][[1]]
      pptid_dpth <- dpths[[pptid]]
      res[i] <- ys[[ptid]][[ptid_dpth - pptid_dpth]]
    }
  }
  res
}

plot(as(tree, 'phylo'))

tree <- randTree(5)
linesdf <- data.frame(x1=rep(NA, tree@nall),
                      x2=NA, y1=NA, y2=NA)
rownames(linesdf) <- tree@all
tree_age <- getAge(tree)
linesdf[tree@all, 'x2'] <- getNdsAge(tree, tree@all,
                                     tree_age=tree_age)
prids <- getNdsSlt(tree, 'prid', tree@all)
linesdf[tree@all, 'x1'] <- linesdf[prids, 'x2']
linesdf[tree@root, c('y1', 'y2')] <- 0
dpths <- sapply(getNdsPtids(tree, tree@all), length) + 1
dpths <- sort(dpths, decreasing=TRUE)
ids <- names(dpths[dpths != 1])
ids <- ids[ids != tree@root]
# start with root, determine possible ys
ptids <- tree@ndlst[[tree@root]][['ptid']]
ys <- calcRootYs(0, ptids)
linesdf[ptids, 'y2'] <- calcY(ys, ptids)
for(id in ids) {
  y2 <- linesdf[id, 'y2']
  ptids <- tree@ndlst[[id]][['ptid']]
  strt <- 1
  for(ptid in ptids) {
    [strt:dpths[[ptid]]]
  }
  cumsum(dpths[ptids])
  ys <- c(ys, calcYs(, ptids))
  linesdf[ptids, 'y2'] <- calcY(ys[ptids], ptids)
}
linesdf[tree@all, 'y1'] <- linesdf[prids, 'y2']
p_data <- data.frame(x=c(linesdf[['x1']],
                         linesdf[['x2']]),
                     y=c(linesdf[['y1']],
                         linesdf[['y2']]),
                     id=rep(rownames(linesdf), 2))
ggplot(p_data, aes(x=x, y=y, group=id)) + geom_line()


if (length (edge) == 2) {
  # run recursively for both
  .cpMkPData (x1, y[1], edge[1], p.env)
  .cpMkPData (x1, y[2], edge[2], p.env)
} else if (length (edge) == 1){
  # find the next node on the edge
  next.node <- p.env$tree$edge[edge,2]
  # get n.children
  n.children <- p.env$tree.stats[[next.node]][['n.children']]
  # get the length of the edge
  x2 <- x1 - p.env$tree$edge.length[edge]
  l.data <- data.frame (x=c(x1, x2), y=c(y, y), edge)
  p.env$p.data <- rbind (p.env$p.data, l.data)
  # run again for next edge
  next.edge <- p.env$tree.stats[[next.node]][['next.edges']]
  if (next.node <= length (p.env$tree$tip.label)) {
    # EDIT THIS TO ADD TIP LABELS
    #t.data <<- data.frame (x=(x2+1/N), y=y,
    #                       label=tree$tip.label[next.node])
  } else {
    # use n.children to determine y and reduce overlap
    y.length <- n.children / p.env$N
    ys <- c (y-(y.length/2), y+(y.length/2))
    # add vertical connector
    l.data <- data.frame (x=c(x2, x2), y=ys, edge)
    p.env$p.data <- rbind (p.env$p.data, l.data)
    # move on to next edge
    .cpMkPData (x2, ys, next.edge, p.env)
  }
} else {
  stop ('Invalid tree')
}