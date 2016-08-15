# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'ndlst-methods\'')
test_that('.getNdPridsFrmLst([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  rnd_tip <- sample(tree@tips, 1)
  res <- treeman:::.getNdPridsFrmLst(ndlst, rnd_tip)
  expect_true(ndlst[[rnd_tip]]['prid'] %in% ids[res])
  expect_true(tree@root %in% ids[res])
  rnd_nd <- sample(tree@nds, 1)
  res <- treeman:::.getNdPridsFrmLst(ndlst, rnd_nd)
  expect_true(ndlst[[rnd_nd]]['prid'] %in% ids[res])
  expect_true(tree@root %in% ids[res])
  res <- treeman:::.getNdPridsFrmLst(ndlst, tree@root)
  expect_true(tree@root == ids[res])
})
test_that('.getNdPrdstFrmLst([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  res <- treeman:::.getNdPrdstFrmLst(ndlst, tree@root)
  expect_that(res, equals(0))
  rnd_nd <- sample(tree@nds[tree@nds != tree@root], 1)
  res <- treeman:::.getNdPrdstFrmLst(ndlst, rnd_nd)
  expect_that(res, is_less_than(tree@age))
})
test_that('.getNdPtidsFrmLst([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  rnd_tip <- sample(tree@tips, 1)
  res <- treeman:::.getNdPtidsFrmLst(ndlst, rnd_tip)
  expect_that(length(res), equals(0))
  intrnls <- tree@nds[tree@nds != tree@root]
  res <- treeman:::.getNdPtidsFrmLst(ndlst, id=tree@root)
  expect_true(all(intrnls %in% ids[res]))
  rnd_nd <- sample(intrnls, 1)
  res <- treeman:::.getNdPtidsFrmLst(ndlst, id=rnd_nd)
  expect_false(all(intrnls %in% ids[res]))
})
test_that('.getNdKidsFrmLst([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  res <- treeman:::.getNdKidsFrmLst(ndlst, id=tree@root)
  expect_true(all(tree@tips %in% ids[res]))
  rnd_nd <- sample(tree@nds[tree@nds != tree@root], 1)
  res <- treeman:::.getNdKidsFrmLst(ndlst, id=rnd_nd)
  expect_false(all(tree@tips %in% ids[res]))
  rnd_tip <- sample(tree@tips, 1)
  res <- treeman:::.getNdKidsFrmLst(ndlst, id=rnd_tip)
  expect_that(length(res), equals(0))
})
test_that('.getNdPDFrmLst([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  res <- treeman:::.getNdPDFrmLst(ndlst, id=tree@root)
  expect_that(res, equals(tree@pd))
  rnd_nd <- sample(tree@nds[tree@nds != tree@root], 1)
  res <- treeman:::.getNdPDFrmLst(ndlst, id=rnd_nd)
  expect_that(res, is_less_than(tree@pd))
})
