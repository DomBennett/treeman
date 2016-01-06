# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'manip-methods\'')
test_that('addTip() works', {
  # random tree + basic stats
  tree <- randTree(10)
  pd_before <- tree@pd
  age_before <- tree@age
  ntips_before <- nTips(tree)
  # add random tip
  sister <- sample(tree@tips, 1)
  sister_age <- getNodeAge(tree, sister)
  parent_age <- getNodeAge(tree, tree@nodelist[[sister]]$prenode)
  start <- runif(min=sister_age, max=parent_age, n=1)
  end <- runif(min=0, max=start, n=1)
  tree <- addTip(tree, id='new_tip', sister=sister, start=start, end=end)
  #viz(tree)
  # test if successful
  expect_that(tree@plytms, is_false())
  expect_that(tree@age, equals(age_before))
  expect_that(nTips(tree), equals(ntips_before + 1))
  expect_that(tree@pd, equals(pd_before + (start-end)))
})