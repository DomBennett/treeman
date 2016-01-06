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
  # add a tip with a specified length.
  sister <- sample(tips, 1)
  sister_age <- getNodeAge(tree, sister)
  parent_age <- getNodeAge(tree, tree@nodelist[[sister]]$prenode)
  start <- runif(min=sister_age, max=parent_age, n=1)
  end <- runif(min=0, max=start, n=1)
  tree <- addTip(tree, id='new_tip', sister=sister, start=start, end=end)
  new_pd <- calcPhyDiv(tree, c(tips, 'new_tip'))
  test_that(new_pd, equals(pd + (start - end)))
})
test_that('calcFairProp() works', {
  #TODO
})
test_that('calcDstMtrx() works', {
  #TODO
})