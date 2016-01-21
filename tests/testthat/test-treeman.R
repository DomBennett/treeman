# LIBS
library(treeman)
library(testthat)

# DATA
tree <- randTree(100)

# RUNNING
context('Testing \'TreeMan Class\'')
test_that('vaildObject() works', {
  res <- validObject(tree)
  expect_that(res, is_true())
  tree@nodelist[['n2']][['id']] <- 'oh oh.... invalid ID'
  expect_error(validObject(tree))
})
test_that('[[ works', {
  node <- sample(names(tree@nodelist), 1)
  node <- tree[[node]]
  expect_that(class(node)[[1]], equals('Node'))
})
test_that('[ works', {
  expect_error(tree['not a valid slot name'])
  expect_that(tree['ntips'], equals(100))
  expect_that(tree['nnodes'], equals(99))
  expect_that(tree['age'], equals(tree@age))
})
test_that('.update() works', {
  tree@nodelist[['t1']][['span']] <- NULL
  tree <- treeman:::.update(tree)
  expect_false(tree@wspn)
})