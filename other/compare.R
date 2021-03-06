# Comparing adding tips to a tree for Phylo and TreeMan

# LIBS
library(ggplot2)

phyloBuilder <- function (start_tree, max_size, sample_size) {
  tree <- ape::read.tree(text=start_tree)
  start <- 1
  repeats <- max_size/sample_size
  timings <- rep(NA, repeats)
  for(i in 1:repeats) {
    timings[i] <- system.time(expr={
      for(j in 1:sample_size) {
        start <- start/2
        id <- paste0('t', length(tree$tip.label)+1)
        nid <- paste0('n', length(tree$tip.label)+2)
        tree <- MoreTreeTools:::addTip(tree, edge=1,
                                       tip.age=0,
                                       node.age=start,
                                       tip.name=id,
                                       node.label=nid)
      }})[3]
  }
  timings
}

treemanBuilder <- function (start_tree, max_size, sample_size) {
  require(treeman)
  tree <- readTree(text=start_tree)
  start <- tree_age <- getAge(tree)
  repeats <- max_size/sample_size
  timings <- rep(NA, repeats)
  for(i in 1:repeats) {
    timings[i] <- system.time(expr={
      for(j in 1:sample_size) {
        id <- paste0('t', tree['ntips']+1)
        tree <- treeman:::addTip(tree=tree, tid=id, sid='t1',
                                 strt_age=start, end=-1,
                                 tree_age=tree_age)
        start <- 1
        tree_age <- tree_age + 1
      }})[3]
  }
  timings
}

# RUN
start_tree <- "(t1:1.0,t2:1.0);"
max_size <- 100  # maximum size of tree
sample_size <- 10
treeman_res <- treemanBuilder(start_tree, max_size, sample_size)
phylo_res <- phyloBuilder(start_tree, max_size, sample_size)

# PLOT
res <- data.frame("Tree_size"=seq(2, max_size, sample_size),
                  "Timing"=c(phylo_res, treeman_res),
                  "Build_method"=rep(c('Phylo', "TreeMan"),
                                     each=max_size/sample_size))
#save(res, file=file.path('other', "compare.Rd"))
p <- ggplot (res, aes(y=Timing, x=Tree_size, colour=Build_method)) + geom_line(size=4)
p + theme_bw() + theme(text=element_text(size=20)) + ylab(label = "Timing(s)") +
  xlab(label = "N. tips")
