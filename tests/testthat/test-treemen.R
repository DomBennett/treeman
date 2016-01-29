# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'TreeMen Class\'')
test_that("TreeMen works", {
  tree <- randTree(10)
  trees <- cTrees(tree, tree, tree)
  expect_true(is(trees, "TreeMen"))
})