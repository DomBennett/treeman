# Comparing adding tips to a tree for Phylo and TreeMan
# D.J. Bennett

# LIBS
library(ggplot2)

addTip <- function(tree, id, sister, start, end,
                   parent_id=paste0("p_", id),
                   tip_taxonym=NULL, parent_taxonym=NULL) {
  updatePre <- function(node) {
    node[['kids']] <- c(node[['kids']], tip[['id']])
    node[['pd']] <- node[['pd']] + tip[['span']]
    node
  }
  tip <- list('id'=id)
  if(!is.null(tip_taxonym)) {
    tip[['taxonym']] <- tip_taxonym
  }
  node <- list('id'=parent_id)
  if(!is.null(parent_taxonym)) {
    node[['taxonym']] <- parent_taxonym
  }
  tip[['span']] <- start - end
  age <- getNodeAge(tree, sister)
  new_sister <- sister <- tree@nodelist[[sister]]
  new_parent <- tree@nodelist[[sister[['prid']]]]
  new_parent[['ptid']] <- new_parent[['ptid']][!new_parent[['ptid']] %in% sister[['id']]]
  new_parent[['ptid']] <- c(new_parent[['ptid']], node[['id']])
  new_sister[['span']] <- start - age
  new_sister[['prid']] <- node[['id']]
  node[['span']] <- sister[['span']] - new_sister[['span']]
  node[['pd']] <- new_sister[['span']] + tip[['span']]
  node[['prdst']] <- sister[['prdst']] - new_sister[['span']]
  node[['prid']] <- sister[['prid']]
  node[['ptid']] <- node[['children']] <- c(tip[['id']], sister[['id']])
  tip[['pd']] <- 0
  tip[['prdst']] <- node[['prdst']] + tip[['span']]
  tip[['prid']] <- node[['id']]
  tree@nodelist[[tip[['id']]]] <- tip
  tree@nodelist[[node[['id']]]] <- node
  tree@nodelist[[new_sister[['id']]]] <- new_sister
  tree@nodelist[[new_parent[['id']]]] <- new_parent
  pres <- getNodePrid(tree, node[['id']])
  tree@nodelist[pres] <- lapply(tree@nodelist[pres],
                                updatePre)
  treeman:::.updateSlots(tree)
}

phyloBuilder <- function (start_tree, max_size, sample_size) {
  require(MoreTreeTools)
  tree <- read.tree(text=start_tree)
  start <- 1
  repeats <- max_size/sample_size
  timings <- rep(NA, repeats)
  for(i in 1:repeats) {
    timings[i] <- system.time(expr={
      for(j in 1:sample_size) {
        start <- start/2
        id <- paste0('t', length(tree$tip.label)+1)
        nid <- paste0('n', length(tree$tip.label)+2)
        tree <- MoreTreeTools:::.addTip__phylo(tree, edge=1,
                                               tip_age=0,
                                               node_age=start,
                                               tip_name=id,
                                               node_name=nid)
      }})[3]
  }
  timings
}

treemanBuilder <- function (start_tree, max_size, sample_size) {
  require(treeman)
  tree <- readTree(text=start_tree)
  start <- tree['age']
  repeats <- max_size/sample_size
  timings <- rep(NA, repeats)
  for(i in 1:repeats) {
    timings[i] <- system.time(expr={
      for(j in 1:sample_size) {
        start <- start/2
        id <- paste0('t', tree['ntips']+1)
        tree <- treeman:::addTip(tree=tree, id=id, sister='t1',
                                 start=start, end=0)
      }})[3]
  }
  timings
}

# RUN
start_tree <- "(t1:1.0,t2:1.0);"
max_size <- 1000  # maximum size of tree
sample_size <- 10
treeman_res <- treemanBuilder(start_tree, max_size, sample_size)
phylo_res <- phyloBuilder(start_tree, max_size, sample_size)

# PLOT
res <- data.frame("Tree_size"=seq(2, max_size, sample_size),
                  "Timing"=c(phylo_res, treeman_res),
                  "Build_method"=rep(c('Phylo', "TreeMan"), each=max_size/sample_size))
#save(res, file=file.path('other', "compare.Rd"))
p <- ggplot (res, aes(y=Timing, x=Tree_size, colour=Build_method)) + geom_line(size=4)
p + theme_bw() + theme(text=element_text(size=20)) + ylab(label = "Timing(s)") +
  xlab(label = "N. tips")
