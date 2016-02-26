# Evolutionary Distinctness Biased Markov Model
# Demonstrating how to simulate a tree with treeman

# LIBS
library(treeman)

# PARAMETER
# balanced tree to start
tree_string <- "((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);"
tree <- readTree(text=tree_string)
iterations <- 200
burnin <- 10
d <- 1
b <- 2  # birth for burnin
b_true <- 1  # birth after burnin
ext <- tree["tips"]

# LOOP
cat('Simulating ....\n')
for(i in 1:iterations) {
  if(length(ext) < 3) {
    stop('Too few tips remaining!')
  }
  if(i > burnin) {
    b <- b_true
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
    tid <- sample(ext, prob=1/fps, size=1)
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
ed <- calcEdgeDiversity(tree_phylo, n.intervals=10)
ed$col <- (log(ed$count) - mean(log(ed$count))) /
  sd(log(ed$count))
p <- chromatophylo(tree_phylo, edge.cols=ed, legend.title='Diversity')
print(p)
