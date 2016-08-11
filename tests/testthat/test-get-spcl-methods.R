# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'get-spcl-methods\'')
test_that('getTreeAge() works', {
  tree <- randTree(10)
  age <- getTreeAge(tree)
  expect_true(age > 0)
  expect_true(age < tree['pd'])
})
test_that('getOtgrp() works', {
  tree <- randTree(10)
  rnd_nd <- sample(tree['nds'][tree['nds'] != tree['root']], 1)
  ingrp <- getNdKids(tree, rnd_nd)
  otgrp <- sample(tree['tips'][!tree['tips'] %in% ingrp], 1)
  res <- getOtgrp(tree, ids=c(ingrp, otgrp))
  expect_that(res, equals(otgrp))
})
test_that("getPrnt() works", {
  tree <- readTree(text="(((A,B),(C,D)),(E,F));")
  prnt <- getPrnt(tree, ids=c("A", "F"))
  expect_true(prnt == tree@root)
})
test_that("getPath() works", {
  tree <- randTree(10)
  pth <- getPath(tree, from="t1", to="t10")
  prnt <- getPrnt(tree, ids=c('t1', "t10"))
  expect_true(prnt %in% pth)
  expect_that(pth[1], equals('t1'))
  expect_that(pth[length(pth)], equals('t10'))
})
test_that("getSubtree() works", {
  tree <- randTree(10)
  subtree <- getSubtree(tree, 'n2')
  expect_that(tree['ntips'], is_more_than(subtree['ntips']))
  expect_that(tree['nnds'], is_more_than(subtree['nnds']))
  expect_that(tree['pd'], is_more_than(subtree['pd']))
  expect_that(tree['age'], is_more_than(subtree['age']))
})