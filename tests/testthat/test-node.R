# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'Node Class\'')
test_that('Node works (rooted + with spns)', {
  tree <- randTree(10)
  nid <- names(tree@nodelist)[sample(1:10, 1)]
  node <- tree[[nid]]
})
test_that('Node works (rooted + w/o spns)', {
  tree <- randTree(10)
  ids <- c( tree['tips'], tree['nodes'])
  #tree <- setNodes(tree, ids=ids, name='span', values=NULL)
  #TODO: remove branch lengths
  nid <- names(tree@nodelist)[sample(1:10, 1)]
  node <- tree[[nid]]
  str(node)
})
test_that('Node works (unrooted  + w/o spns))', {
  tree <- randTree(10)
  #TODO: remove branch lengths
  #TODO: remove root
  nid <- names(tree@nodelist)[sample(1:10, 1)]
  node <- tree[[nid]]
})