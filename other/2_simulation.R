# Evolutionary Distinctness Biased Markov Model
# Demonstrating how to simulate a tree with treeman

# LIBS
library(treeman)

# PARAMETER
# balanced tree to start
tree_string <- "((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);"
tree_age <- 2
tree <- readTree(text=tree_string)
iterations <- 200
burnin <- 10
d <- 1
b <- 3  # birth for burnin
b_true <- 1  # birth after burnin
extnt <- tree["tips"]
extnct <- NULL

# LOOP
cat('Simulating ....\n')
for(i in 1:iterations) {
  if(length(extnt) < 3) {
    stop('Too few extant tips remaining!\nRun again or increase birth before burnin.')
  }
  if(i > burnin) {
    b <- b_true
  }
  cat('.... i=[', i, ']\n', sep='')
  # calculate fair proportion
  fps <- calcPrtFrPrp(tree, tids=extnt, ignr=extnct)
  # add/remove based on b and d
  to_add <- sample(c(TRUE, FALSE), size=1, prob=c(b,d))
  if(to_add) {
    sid <- sample(extnt, prob=1/fps, size=1)  # sister ID of the new tip
    tid <- paste0('t', i)  # new tip ID
    tree <- addTip(tree, tid=tid, sid=sid, strt_age=0, end_age=0,
                   tree_age=tree_age)
    extnt <- c(extnt, tid)
  } else {
    tid <- sample(extnt, prob=1/fps, size=1)
    extnct <- c(extnct, tid)
    extnt <- extnt[extnt != tid]
  }
  # grow tree
  spns <- getNdsSlt(tree, slt_nm="spn", ids=extnt)
  tree <- setNdsSpn(tree, ids=extnt, vals=spns+1)
  tree_age <- tree_age + 1
}
cat('Done.\n')

# VIZ
# plot the fossil record tree, and the reconstructed phylogeny
extnt_tree <- rmTips(tree, tids=getDcsd(tree))
par(mfrow=c(1,2))
plot(as(tree, 'phylo'), show.tip.label=FALSE, main='Fossil record')
plot(as(extnt_tree, 'phylo'), main='Reconstructed')
