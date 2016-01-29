# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'calc-methods\'')
test_that('calcDstTrp() works', {
  tree_1 <- randTree(10)
  tree_2 <- randTree(10)
  res <- calcDstTrp(tree_1, tree_2, nrmlsd=TRUE)
  expect_that(res, is_more_than(0))
  res <- calcDstTrp(tree_1, tree_1, nrmlsd=TRUE)
  expect_that(res, equals(0))
})
test_that('calcOvrlp() works', {
  tree <- readTree(tree_string="((t1:1.0,t2:1.0):1.0,t3:1.0);")
  ovrlp <- calcOvrlp(tree, ids_1=tree['tips'], ids_2=c('t3'), nrmlsd=TRUE)
  expect_that(ovrlp, equals(1/4))
})
test_that('calcDstBLD() works', {
  tree_1 <- readTree(tree_string="((t1:1.0,t2:1.0):1.0,t3:1.0);")
  tree_2 <- readTree(tree_string="((t3:1.0,t2:1.0):1.0,t1:1.0);")
  d <- calcDstBLD(tree_1, tree_2, TRUE)
  expect_that(d, equals(1))
  d <- calcDstBLD(tree_1, tree_1, TRUE)
  expect_that(d, equals(0))
})
test_that('calcDstRF() works', {
  tree_1 <- readTree(tree_string="((t1,t2),t3);")
  tree_2 <- readTree(tree_string="((t3,t2),t1);")
  d <- calcDstRF(tree_1, tree_2, TRUE)
  expect_that(d, equals(1))
  d <- calcDstRF(tree_1, tree_1, TRUE)
  expect_that(d, equals(0))
})
test_that('calcPhyDv() works', {
  tree <- randTree(10)
  tips <- sample(tree['tips'], 3)
  pd <- calcPhyDv(tree, tips)
  parent <- getParent(tree, ids=tips)
  test_that(pd, is_less_than(tree@nodelist[[parent]][['pd']]))
  # add a tip with a specified length.
  sister <- sample(tips, 1)
  sister_age <- getNodeAge(tree, sister)
  parent_age <- getNodeAge(tree, tree@nodelist[[sister]][['prid']])
  start <- runif(min=sister_age, max=parent_age, n=1)
  end <- runif(min=0, max=start, n=1)
  tree <- addTip(tree, id='new_tip', sister=sister, start=start, end=end)
  new_pd <- calcPhyDv(tree, c(tips, 'new_tip'))
  test_that(new_pd, equals(pd + (start - end)))
})
test_that('calcFrPrp() works', {
  tree <- randTree(10)
  ed_values <- calcFrPrp(tree, tree['tips'])
  expect_that(sum(ed_values), equals(tree['pd']))
})
test_that('calcDstMtrx() works', {
  tree <- randTree(10)
  ids <- tree['all']
  dmtrx <- calcDstMtrx(tree, ids)
  rndnd <- sample(ids, 1)
  expect_that(dmtrx[rndnd, rndnd], equals(0))
  expect_that(sum(dmtrx['n1', ] == tree['age']), is_more_than(0))
})