# Evolutionary Distinctness Biased Markov Model
# Demonstrating how to simulate a tree with treeman
# N.B. simulation may fail if all tips go extinct

# LIBS
library(treeman)

# PARAMETER
# balanced tree to start
tree_string <- "((t1:1.0,t2:1.0)1:0,(t3:1.0,t4:1.0):1.0);"
tree <- readTree(text=tree_string)
iterations <- 10
b <- 2
d <- 1
exc <- NULL
ext <- tree["tips"]

# LOOP
for(i in 1:iterations) {
  # calculate fair proportion
  fps <- calcFrPrp(tree, ext)
  # add/remove based on b and d
  to_add <- sample(c(TRUE, FALSE), size=1, prob=c(b,d))
  if(to_add) {
    sister <- sample(ext, prob=fps, size=1)
    id <- paste0('t', tree['ntips']+1)
    tree <- addTip(tree, id=id, sister=sister, start=0, end=0)
  } else {
    exc <- c(exc, sample(ext, prob=1/fps, size=1))
  }
  # grow tree
  ext <- tree['tips'][!tree['tips'] %in% exc]
  spans <- getNodesSlot(tree, name="span", ids=ext)
  tree <- setNodesSpan(tree, ids=ext, vals=spans+1)
}

# VIZ
viz(tree)
