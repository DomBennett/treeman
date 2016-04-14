# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'set-methods\'')
test_that('set_ID() works', {
  tree <- randTree(100)
  vals <- paste0('new_id_', 1:100)
  ids <- tree['tips']
  tree <- setNdsID(tree, ids=ids, vals=vals)
  expect_true(all(tree['tips'] == vals))
  expect_true(all(tree[['n1']]['kids'] %in% vals))
  tree <- setNdID(tree, id='new_id_1', val='t1')
  expect_error(tree[['new_id_1']])
})
test_that('setNdSpn() works', {
  tree <- randTree(10)
  before_age <- tree['age']
  before_pd <- tree['pd']
  ids <- tree['all'][tree['all'] != tree['root']]
  id <- sample(ids, 1)
  before_prdst <- tree[[id]]['prdst']
  val <- tree[[id]]['spn']/2
  tree <- setNdSpn(tree, id=id, val=val)
  expect_that(tree['pd'] + val, equals(before_pd))
  expect_that(tree[[id]]['prdst'] + val, equals(before_prdst))
})
test_that('setNdsSpn() works', {
  tree <- randTree(10)
  before_pd <- tree['pd']
  before_age <- tree['age']
  ids <- tree['all'][tree['all'] != tree['root']]
  vals <- getNdsSlt(tree, name='spn', ids=ids)
  vals <- vals/2
  tree <- setNdsSpn(tree, ids=ids, vals=vals)
  expect_that(tree['pd']*2, equals(before_pd))
  expect_that(tree['age']*2, equals(before_age))
  tree <- setNdsSpn(tree, ids=ids, vals=NULL)
  expect_false(tree['wspn'])
})
test_that('setPD() works', {
  tree <- randTree(10)
  tree <- setPD(tree, val=1)
  expect_that(tree['pd'], equals(1))
})
test_that('setAge() works', {
  tree <- randTree(10)
  tree <- setAge(tree, val=1)
  expect_that(tree['age'], equals(1))
})
test_that('setTol() works', {
  tree <- randTree(10)
  before <- length(tree@ext)
  tree <- setTol(tree, tol=tree['age'])
  expect_that(before, is_less_than(length(tree@ext)))
})