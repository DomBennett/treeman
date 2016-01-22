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
test_that('set_Span() works', {
  before <- tree['pd']
  val <- tree[['t1']]['span']/2
  tree <- setNodeSpan(tree, id='t1', val=val)
  expect_that(tree['pd'] + val, equals(before))
  before <- tree['pd']
  ids <- c(tree['tips'], tree['nodes'])
  vals <- getNodesSlot(tree, name='span', ids=ids)
  vals <- vals/2
  tree <- setNodesSpan(tree, ids=ids, vals=vals)
  expect_that(tree['pd']*2, equals(before))
  tree <- setNodesSpan(tree, ids=ids, vals=NULL)
  expect_false(tree['wspn'])
})
test_that('setPD() works', {
  tree <- setPD(tree, val=1)
  expect_that(tree['pd'], equals(1))
})
test_that('setAge() works', {
  tree <- setAge(tree, val=1)
  expect_that(tree['age'], equals(1))
})
test_that('setTol() works', {
  before <- length(tree@ext)
  tree <- setTol(tree, tol=tree['age'])
  expect_that(before, is_less_than(length(tree@ext)))
})