# TODO:
# -- read/write newick
library(treeman)

tree <- readTree("~/Desktop/bininda.tre")


tips(tree)
age(tree)
calcFairProp(tree, sample(tips(tree), 10))
calcDstMtrx(tree)

tree[['n1']]
tree[['n5']]
tree[['n4']]

trstr <- "((A:1.0,B:1.0):1.0,(C:1.0,(E:1.0,D:1.0):1.0):1.0);"
tree <- readTree(tree_string=trstr)
viz(tree)
