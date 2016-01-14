# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'read-write-methods\'')
test_that('readTree() works', {
  tree_string <- "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5);"
  tree <- readTree(tree_string=tree_string)
})