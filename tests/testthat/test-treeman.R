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
test_that('nTips() works', {
  expect_that(length(tree@tips), equals(nTips(tree)))
})
test_that('nodes() works', {
  nodes <- nodes(tree)
  expect_that(length(nodes), equals(99))
})
test_that('nNodes() works', {
  expect_that(length(tree@nodes), equals(nNodes(tree)))
})
test_that('age() works', {
  expect_that(tree@age, equals(age(tree)))
})
test_that('pd() works', {
  expect_that(tree@pd, equals(pd(tree)))
})
test_that('extant() works', {
  expect_that(tree@extant, equals(extant(tree)))
})
test_that('extinct() works', {
  expect_that(tree@extinct, equals(extinct(tree)))
})
test_that('rootNode() works', {
  expect_that(tree@root, equals(rootNode(tree)))
})
test_that('ultrmtrc() works', {
  expect_that(tree@ultrmtrc, equals(ultrmtrc(tree)))
})
test_that('plytms() works', {
  expect_that(tree@plytms, equals(plytms(tree)))
})