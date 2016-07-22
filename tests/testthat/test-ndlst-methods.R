# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'ndlst-methods\'')
test_that('.getNdPrids([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  rnd_tip <- sample(tree@tips, 1)
  res <- treeman:::.getNdPrids(ndlst, rnd_tip)
  expect_true(ndlst[[rnd_tip]]['prid'] %in% ids[res])
  expect_true(tree@root %in% ids[res])
  rnd_nd <- sample(tree@nds, 1)
  res <- treeman:::.getNdPrids(ndlst, rnd_nd)
  expect_true(ndlst[[rnd_nd]]['prid'] %in% ids[res])
  expect_true(tree@root %in% ids[res])
  res <- treeman:::.getNdPrids(ndlst, tree@root)
  expect_true(tree@root == ids[res])
})
test_that('.getNdPrdst([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  res <- treeman:::.getNdPrdst(ndlst, tree@root)
  expect_that(res, equals(0))
  rnd_nd <- sample(tree@nds[tree@nds != tree@root], 1)
  res <- treeman:::.getNdPrdst(ndlst, rnd_nd)
  expect_that(res, is_less_than(tree@age))
})
test_that('.getNdPtids([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  rnd_tip <- sample(tree@tips, 1)
  res <- treeman:::.getNdPtids(ndlst, rnd_tip)
  expect_that(length(res), equals(0))
  intrnls <- tree@nds[tree@nds != tree@root]
  res <- treeman:::.getNdPtids(ndlst, id=tree@root)
  expect_true(all(intrnls %in% ids[res]))
  rnd_nd <- sample(intrnls, 1)
  res <- treeman:::.getNdPtids(ndlst, id=rnd_nd)
  expect_false(all(intrnls %in% ids[res]))
})
test_that('.getNdKids([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  res <- treeman:::.getNdKids(ndlst, id=tree@root)
  expect_true(all(tree@tips %in% ids[res]))
  rnd_nd <- sample(tree@nds[tree@nds != tree@root], 1)
  res <- treeman:::.getNdKids(ndlst, id=rnd_nd)
  expect_false(all(tree@tips %in% ids[res]))
  rnd_tip <- sample(tree@tips, 1)
  res <- treeman:::.getNdKids(ndlst, id=rnd_tip)
  expect_that(length(res), equals(0))
})
test_that('.getNdPD([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  res <- treeman:::.getNdPD(ndlst, id=tree@root)
  expect_that(res, equals(tree@pd))
  rnd_nd <- sample(tree@nds[tree@nds != tree@root], 1)
  res <- treeman:::.getNdPD(ndlst, id=rnd_nd)
  expect_that(res, is_less_than(tree@pd))
})
test_that('.getTreeAge([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  # tests
  res <- treeman:::.getTreeAge(ndlst)
  expect_that(res, equals(tree@age))
})
