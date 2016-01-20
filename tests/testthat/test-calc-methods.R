# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'calc-methods\'')
test_that('calcPhyDv() works', {
  tree <- randTree(10)
  tips <- sample(tree@tips, 3)
  pd <- calcPhyDv(tree, tips)
  parent <- getParent(tree, nodes=tips)
  test_that(pd, is_less_than(tree@nodelist[[parent]]$pd))
  # add a tip with a specified length.
  sister <- sample(tips, 1)
  sister_age <- getNodeAge(tree, sister)
  parent_age <- getNodeAge(tree, tree@nodelist[[sister]]$prenode)
  start <- runif(min=sister_age, max=parent_age, n=1)
  end <- runif(min=0, max=start, n=1)
  tree <- addTip(tree, id='new_tip', sister=sister, start=start, end=end)
  new_pd <- calcPhyDv(tree, c(tips, 'new_tip'))
  test_that(new_pd, equals(pd + (start - end)))
})
test_that('calcFrPrp() works', {
  tree <- randTree(10)
  ed_values <- calcFrPrp(tree, tips(tree))
  expect_that(sum(ed_values), equals(pd(tree)))
})
test_that('calcDstMtrx() works', {
  tree <- randTree(10)
  dmtrx <- calcDstMtrx(tree)
  rndnd <- sample(c(nodes(tree), tips(tree)), 1)
  expect_that(dmtrx[rndnd, rndnd], equals(0))
  expect_that(sum(dmtrx['n1', ] == age(tree)), is_more_than(0))
})