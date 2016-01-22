# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'Node Class\'')
test_that('.newNode works', {
  node <- treeman:::.newNode(randTree(10), 'n1')
  expect_that(class(node)[[1]], equals('Node'))
})
test_that('Node works (rooted + with spns)', {
  tree <- randTree(10)
  nid <- names(tree@nodelist)[sample(1:10, 1)]
  node <- tree[[nid]]
})
test_that('Node works (rooted + w/o spns)', {
  tree <- randTree(10)
  tree <- setNodesSpan(tree, ids=NULL, vals=NULL)
  nid <- names(tree@nodelist)[sample(1:10, 1)]
  node <- tree[[nid]]
})
test_that('Node works (unrooted  + w/o spns))', {
  tree <- randTree(10)
  #TODO: remove branch lengths
  #TODO: remove root
  nid <- names(tree@nodelist)[sample(1:10, 1)]
  node <- tree[[nid]]
})