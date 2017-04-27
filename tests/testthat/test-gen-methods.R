# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'gen-methods\'')
test_that('randTree() works', {
  ns <- sample(3:100, 5)
  for(n in ns) {
    tree <- randTree(n)
    expect_that(tree['ntips'], equals(n))
    tree_age <- getAge(tree)
    expect_that(tree['pd'], is_more_than(tree_age))
  }
})
test_that('blncdTree() works', {
  ns <- sample(3:100, 5)
  for(n in ns) {
    tree <- blncdTree(n)
    expect_that(tree['ntips'], equals(n))
    tree_age <- getAge(tree)
    expect_that(tree['pd'], is_more_than(tree_age))
  }
})
test_that('unblncdTree() works', {
  ns <- sample(3:100, 5)
  for(n in ns) {
    tree <- unblncdTree(n)
    expect_that(tree['ntips'], equals(n))
    tree_age <- getAge(tree)
    expect_that(tree['pd'], is_more_than(tree_age))
  }
})