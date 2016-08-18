# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'cnvrt-methods\'')
test_that('TreeMan-to-phylo works', {
  tree <- randTree(100)
  tree <- as(tree, "phylo")
  expect_equal(class(tree)[1], "phylo")
})
test_that('phylo-to-TreeMan works', {
  tree <- ape::rtree(100)
  tree <- as(tree, "TreeMan")
  expect_equal(class(tree)[1], "TreeMan")
})