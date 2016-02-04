# Comparing adding tips to a tree for Phylo and TreeMan
# D.J. Bennett

# LIBS
library(ggplot2)

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
max_size <- 500  # maximum size of tree
sample_size <- 10
phylo_res <- phyloBuilder(start_tree, max_size, sample_size)
treeman_res <- treemanBuilder(start_tree, max_size, sample_size)

# PLOT
res <- data.frame("Tree_size"=seq(2, max_size, sample_size),
                  "Timing"=c(phylo_res, treeman_res),
                  "Build_method"=rep(c('Phylo', "TreeMan"), each=max_size/sample_size))
save(res, file=file.path(getwd(), "compare_results.Rd"))
p <- ggplot (res, aes(y=Timing, x=Tree_size, colour=Build_method)) + geom_line(size=4)
p + theme_bw() + theme(text=element_text(size=20))
