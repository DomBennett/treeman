# Evolutionary Distinctness Biased Markov Model
# Demonstrating how to simulate a tree with treeman

# LIBS
library(treeman)

# PARAMETER
# balanced tree to start
tree_string <- "((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);"
tree <- readTree(text=tree_string)
iterations <- 1000
burnin <- iterations*.25
b <- 1
d <- 0  # death for burnin
d_true <- 1  # death after burnin
ext <- tree["tips"]

# LOOP
cat('Simulating ....\n')
for(i in 1:iterations) {
  if(length(ext) < 3) {
    stop('Too few tips remaining!')
  }
  if(i > burnin) {
    d <- d_true
  }
  cat('.... i=[', i, ']\n', sep='')
  # calculate fair proportion
  fps <- calcFrPrp(tree, ext)
  # add/remove based on b and d
  to_add <- sample(c(TRUE, FALSE), size=1, prob=c(b,d))
  if(to_add) {
    sid <- sample(ext, prob=1/fps, size=1)  # sister ID of the new tip
    tid <- paste0('t', i)  # new tip ID
    tree <- treeman::addTip(tree, tid=tid, sid=sid, start=0, end=0)
  } else {
    tid <- sample(ext, prob=fps, size=1)
    tree <- rmTip(tree, tid=tid)
  }
  # grow tree
  ext <- tree['tips']
  spans <- getNodesSlot(tree, name="span", ids=ext)
  tree <- setNodesSpan(tree, ids=ext, vals=spans+1)
}
cat('Done.\n')

# VIZ
library(MoreTreeTools)
tree_phylo <- as(tree, 'phylo')
# plot simulated tree with edges coloured by proximate diversity
tree_phylo$edge.label <- paste0 ('edge_', 1:nrow(tree_phylo$edge))
# intervals are used to calculate the colour of the branch
# diversity is the number of descedents of a branch with an interval
# here we ensure the interval has 20 timesteps
ed <- calcEdgeDiversity(tree_phylo, n.intervals=round(iterations/20))
ed$col <- (log(ed$count) - mean(log(ed$count))) /
  sd(log(ed$count))
chromatophylo(tree_phylo, edge.cols=ed, legend.title='Diversity')
