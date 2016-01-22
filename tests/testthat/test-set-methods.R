# LIBS
library(treeman)
library(testthat)

# DATA
tree <- randTree(100)

# RUNNING
context('Testing \'set-methods\'')
test_that('set_ID() works', {
  vals <- paste0('new_id_', 1:100)
  ids <- tree['tips']
  tree <- setNodesID(tree, ids=ids, vals=vals)
  expect_true(all(tree['tips'] == vals))
  expect_true(all(tree[['n1']]['children'] %in% vals))
  tree <- setNodeID(tree, id='new_id_1', val='t1')
  expect_error(tree[['new_id_1']])
})
test_that('setTol() works', {
  before <- length(tree@ext)
  tree <- setTol(tree, tol=tree['age'])
  expect_that(before, is_less_than(length(tree@ext)))
})