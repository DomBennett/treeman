# LIBS
library(treeman)
library(testthat)

# DATA
tree <- randTree(100)

# RUNNING
context('Testing \'TreeMan Class\'')
test_that('.updateSlots() works', {
  tree@nodelist[['t1']][['span']] <- NULL
  tree <- treeman:::.updateSlots(tree)
  expect_false(tree@wspn)
})