# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'calc-methods\'')
test_that('calc_Blnc() works', {
  tree <- readTree(text="((t1:1.0,t2:1.0):1.0,(t3:1.0,t4:1.0):1.0);")
  bstats <- calcNdsBlnc(tree, tree['nds'])
  expect_that(sum(bstats), equals(0))
})
test_that('calcDstTrp() works', {
  tree_1 <- randTree(10)
  tree_2 <- setNdsID(tree_1, tree_1['tips'],
                     sample(tree_1['tips']))
  res <- calcDstTrp(tree_1, tree_2, nrmlsd=TRUE)
  expect_that(res, is_more_than(0))
  res <- calcDstTrp(tree_1, tree_1, nrmlsd=TRUE)
  expect_that(res, equals(0))
})
test_that('calcOvrlp() works', {
  tree <- readTree(text="((t1:1.0,t2:1.0):1.0,t3:1.0);")
  ovrlp <- calcOvrlp(tree, ids_1=tree['tips'], ids_2=c('t3'), nrmlsd=TRUE)
  expect_that(ovrlp, equals(1/4))
})
test_that('calcDstBLD() works', {
  tree_1 <- readTree(text="((t1:1.0,t2:1.0):1.0,t3:1.0);")
  tree_2 <- readTree(text="((t3:1.0,t2:1.0):1.0,t1:1.0);")
  d <- calcDstBLD(tree_1, tree_2, TRUE)
  expect_that(d, equals(1))
  d <- calcDstBLD(tree_1, tree_1, TRUE)
  expect_that(d, equals(0))
})
test_that('calcDstRF() works', {
  tree_1 <- readTree(text="((t1,t2),t3);")
  tree_2 <- readTree(text="((t3,t2),t1);")
  d <- calcDstRF(tree_1, tree_2, nrmlsd=TRUE)
  expect_that(d, equals(1))
  d <- calcDstRF(tree_1, tree_1, TRUE)
  expect_that(d, equals(0))
})
test_that('calcPhyDv() works', {
  tree <- randTree(10)
  tips <- sample(tree['tips'], 3)
  pd <- calcPhyDv(tree, tips)
  parent <- getPrnt(tree, ids=tips)
  test_that(pd, is_less_than(tree@ndlst[[parent]][['pd']]))
})
test_that('calcFrPrp() works', {
  tree <- randTree(10)
  sum(sapply(tree@ndlst, function(x) x[['spn']]))
  ed_values <- calcFrPrp(tree, tree['tips'])
  expect_that(sum(ed_values), equals(tree['pd']))
})
test_that('calcDstMtrx() works', {
  tree <- randTree(10)
  ids <- tree['all']
  dmtrx <- calcDstMtrx(tree, ids)
  rndnd <- sample(ids, 1)
  expect_that(dmtrx[rndnd, rndnd], equals(0))
  tree_age <- getAge(tree)
  expect_that(sum(dmtrx['n1', ] == tree_age), is_more_than(0))
})