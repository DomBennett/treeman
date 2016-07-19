# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'check-methods\'')
test_that("checkTreeMen() works", {
  tree <- randTree(10)
  trees <- cTrees(tree, tree, tree)
  trees@treelist[[4]] <- "this is not a tree"
  res <- checkTreeMen(trees)
  expect_false(res)
})
test_that("checkTreeMan() works", {
  # test missing node
  tree <- randTree(10)
  tree@ndlst <- tree@ndlst[-(sample(1:10, 1))]
  nonode <- checkTreeMan(tree)
  expect_false(nonode)
  # test absent prid
  tree <- randTree(10)
  rnd <- sample(tree@all, 1)
  tree@ndlst[[rnd]][['prid']] <- NULL
  absent_prid <- checkTreeMan(tree)
  expect_false(absent_prid)
  # test missing root
  tree <- randTree(10)
  tree@root <- 'x'
  noroot <- checkTreeMan(tree)
  expect_false(noroot)
  # test self-reference
  tree <- randTree(10)
  rnd <- sample(tree@nds[tree@nds != tree@root], 1)
  tree@ndlst[[rnd]][['ptid']][1] <- 'x'
  tree@ndlst[[rnd]][['prid']] <- 'x'
  self_reference <- checkTreeMan(tree)
  expect_false(self_reference)
  # test circularity
  rnd <- sample(tree@nds[tree@nds != tree@root], 1)
  tree@ndlst[[rnd]][['ptid']][1] <- tree@root
  circularity <- checkTreeMan(tree)
  expect_false(circularity)
})