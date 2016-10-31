# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'update-methods\'')
test_that('updateTree() works', {
  tree <- randTree(100)
  tree@ndlst[['t1']][['spn']] <- NULL
  expect_error(updateTree(tree))
  tree@ndlst[['t1']][['spn']] <- 1000
  tree <- updateTree(tree)
  expect_that(tree['age'], equals(1000))
  expect_that(tree['pd'], is_more_than(1000))
})
test_that('downdateTree() works', {
  tree <- randTree(100)
  tree <- downdateTree(tree)
  expect_false(tree['updtd'])
  # ndmtrx is null
  expect_true(is.null(tree@ndmtrx))
})