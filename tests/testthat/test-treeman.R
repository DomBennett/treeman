# LIBS
library(treeman)
library(testthat)

# DATA
tree <- randTree(100, wndmtrx=sample(c(TRUE, FALSE), 1))

# RUNNING
context('Testing \'TreeMan Class\'')
test_that('vaildObject() works', {
  expect_true(validObject(tree))
  tree@ndlst[['n2']][['id']] <- 'oh oh.... invalid ID'
  expect_error(validObject(tree))
})
test_that('[[ works', {
  nd <- sample(names(tree@ndlst), 1)
  nd <- tree[[nd]]
  expect_that(class(nd)[[1]], equals('Node'))
})
test_that('[ works', {
  expect_error(tree['not a valid slot name'])
  expect_type(tree['age'], 'double')
  expect_type(tree['ultr'], 'logical')
  expect_that(tree['ntips'], equals(100))
  expect_that(tree['nnds'], equals(99))
})