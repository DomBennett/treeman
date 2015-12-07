# LIBS
library(treeman)
library(testthat)

# DATA
tree <- randTree(100)

# RUNNING
context('Testing \'TreeMan-methods\'')
test_that('vaildObject() works', {
  res <- validObject(tree)
  expect_that(res, is_true())
  tree@nodelist[['n2']]$id <- 'oh oh.... invalid ID'
  expect_error(validObject(tree))
})
test_that('[[ works', {
  node <- sample(names(tree@nodelist), 1)
  node <- tree[[node]]
  expect_that(class(node)[[1]], equals('Node'))
})
test_that('tips() works', {
  tips <- tips(tree)
  expect_that(length(tips), equals(100))
})
test_that('tips<- works', {
  tips(tree)[1] <- 'a new name'
  new_tip <- tree[['a new name']]
  expect_that(class(new_tip)[[1]], equals('Node'))
})
test_that('nTips() works', {
  expect_that(length(tree@tips), equals(nTips(tree)))
})
test_that('nodes() works', {
  nodes <- nodes(tree)
  expect_that(length(nodes), equals(99))
})
test_that('nodes<- works', {
  nodes(tree)[1] <- 'a new name'
  new_node <- tree[['a new name']]
  expect_that(class(new_node)[[1]], equals('Node'))
})
test_that('nNodes() works', {
  expect_that(length(tree@nodes), equals(nNodes(tree)))
})
#TODO: age, pd, extant, extinct, rootNode, ultrmtrc, plytms, setTol