# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'ndlst-methods\'')
test_that('.getPrids([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  rnd_tip <- sample(tree@tips, 1)
  res <- treeman:::.getPrids(ndlst, rnd_tip)
  expect_true(ndlst[[rnd_tip]]['prid'] %in% ids[res])
  expect_true(tree@root %in% ids[res])
  rnd_nd <- sample(tree@nds, 1)
  res <- treeman:::.getPrids(ndlst, rnd_nd)
  expect_true(ndlst[[rnd_nd]]['prid'] %in% ids[res])
  expect_true(tree@root %in% ids[res])
  res <- treeman:::.getPrids(ndlst, tree@root)
  expect_true(tree@root == ids[res])
})
test_that('.getPrdst([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  res <- treeman:::.getPrdst(ndlst, tree@root)
  expect_that(res, equals(0))
  rnd_nd <- sample(tree@nds[tree@nds != tree@root], 1)
  res <- treeman:::.getPrdst(ndlst, rnd_nd)
  expect_that(res, is_less_than(tree@age))
})
test_that('.getPtids([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  rnd_tip <- sample(tree@tips, 1)
  res <- treeman:::.getPtids(ndlst, rnd_tip)
  expect_that(length(res), equals(0))
  intrnls <- tree@nds[tree@nds != tree@root]
  res <- treeman:::.getPtids(ndlst, id=tree@root)
  expect_true(all(intrnls %in% ids[res]))
  rnd_nd <- sample(intrnls, 1)
  res <- treeman:::.getPtids(ndlst, id=rnd_nd)
  expect_false(all(intrnls %in% ids[res]))
})
test_that('.getKids([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  res <- treeman:::.getKids(ndlst, id=tree@root)
  expect_true(all(tree@tips %in% ids[res]))
  rnd_nd <- sample(tree@nds[tree@nds != tree@root], 1)
  res <- treeman:::.getKids(ndlst, id=rnd_nd)
  expect_false(all(tree@tips %in% ids[res]))
  rnd_tip <- sample(tree@tips, 1)
  res <- treeman:::.getKids(ndlst, id=rnd_tip)
  expect_that(length(res), equals(0))
})
test_that('.getPD([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  res <- treeman:::.getPD(ndlst, id=tree@root)
  expect_that(res, equals(tree@pd))
  rnd_nd <- sample(tree@nds[tree@nds != tree@root], 1)
  res <- treeman:::.getPD(ndlst, id=rnd_nd)
  expect_that(res, is_less_than(tree@pd))
})
test_that('.getAge([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # tests
  res <- treeman:::.getAge(ndlst)
  expect_true(res == tree@age)
})