# LIBS
library(treeman)
library(testthat)

# DATA
tree <- randTree(100)

# RUNNING
context('Testing \'set-methods\'')
# test_that('nodes<- works', {
#   nodes(tree)[1] <- 'a new name'
#   new_node <- tree[['a new name']]
#   expect_that(class(new_node)[[1]], equals('Node'))
# })
# test_that('tips<- works', {
#   tips(tree)[1] <- 'a new name'
#   new_tip <- tree[['a new name']]
#   expect_that(class(new_tip)[[1]], equals('Node'))
# })