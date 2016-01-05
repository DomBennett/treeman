# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'calc-methods\'')
test_that('calcPhyDiv() works', {
  tree <- randTree(10)
  tips <- sample(tree@tips, 3)
  pd <- calcPhyDiv(tree, tips)
  parent <- getParent(tree, nodes=tips)
  test_that(pd, is_less_than(tree@nodelist[[parent]]$pd))
  #TODO use addTip to add a tip with a specified length.
})
test_that('calcFairProp() works', {
  #TODO
})
test_that('calcDstMtrx() works', {
  #TODO
})