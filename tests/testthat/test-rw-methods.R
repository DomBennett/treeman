# LIBS
library(treeman)
library(testthat)

# DATA
data(mammals)

# RUNNING
context('Testing \'read-write-methods\'')
test_that('readTree([w/ spans]) works', {
  text <- "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5);"
  tree <- readTree(text=text)
  expect_that(tree['pd'], equals(0.1+0.2+0.3+0.4+0.5))
  expect_that(tree['ntips'], equals(4))
  expect_that(tree[['C']]['prdst'], equals(0.3+0.5))
})
test_that('readTree([w/o spans]) works', {
  text <- "(A,B,(C,D)E);"
  tree <- readTree(text=text)
  expect_that(tree['ntips'], equals(4))
  expect_false(tree['wspn'])
})
test_that('writeTree() works', {
  t1 <- randTree(100)
  writeTree(t1, 'test.tre')
  t2 <- readTree('test.tre')
  expect_that(t1['ntips'], equals(t2['ntips']))
  expect_that(t1['nnds'], equals(t2['nnds']))
  expect_that(t1['pd'], equals(t2['pd']))
  t1_age <- getAge(t1)
  t2_age <- getAge(t2)
  expect_that(t1_age, equals(t2_age))
})
test_that('writeTrmn([for TreeMan]) works', {
  t1 <- randTree(100)
  writeTrmn(t1, 'test.trmn')
  t2 <- readTrmn('test.trmn')
  expect_that(t1['ntips'], equals(t2['ntips']))
  expect_that(t1['nnds'], equals(t2['nnds']))
  expect_that(t1['pd'], equals(t2['pd']))
  t1_age <- getAge(t1)
  t2_age <- getAge(t2)
  expect_that(t1_age, equals(t2_age))
  # test example with taxonomy
  data(mammals)
  ape_id <- getPrnt(mammals, ids=c('Homo_sapiens', 'Hylobates_concolor'))
  tree <- getSubtree(mammals, id=ape_id)
  writeTrmn(tree, file='test.trmn')
  tree <- readTrmn('test.trmn')
  expect_true(tree@wtxnyms)
})
test_that('writeTrmn([for TreeMen]) works', {
  t1 <- randTree(100)
  ape_id <- getPrnt(mammals, ids=c('Homo_sapiens', 'Hylobates_concolor'))
  t2 <- getSubtree(mammals, id=ape_id)
  t3 <- randTree(100)
  t3 <- setNdsSpn(t3, t3['all'], vals=0)
  trs <- cTrees(t1, t2, t3)
  writeTrmn(trs, file='test.trmn')
  expect_error(writeTrmn('not_a_tree', file='test.trmn'))
})
test_that('readTrmn([for TreeMan]) works', {
  tree <- randTree(100)
  writeTrmn(tree, 'test.trmn')
  t1 <- readTrmn('test.trmn')
  tree@wspn <- FALSE
  writeTrmn(tree, 'test.trmn')
  t2 <- readTrmn('test.trmn')
  expect_true(t1@wspn)
  expect_true(!t2@wspn)
})
test_that('readTrmn([for TreeMen]) works', {
  t1 <- randTree(100)
  ape_id <- getPrnt(mammals, ids=c('Homo_sapiens', 'Hylobates_concolor'))
  t2 <- getSubtree(mammals, id=ape_id)
  t3 <- randTree(100)
  t3 <- setNdsSpn(t3, t3['all'], vals=0)
  trs <- cTrees(t1, t2, t3)
  writeTrmn(trs, file='test.trmn')
  remove(trs)
  trs <- readTrmn(file='test.trmn')
  expect_true(is(trs) == 'TreeMen')
  expect_true(trs[[1]]['ntips'] == 100)
  expect_true(trs[[1]]['wspn'])
  expect_false(trs[[1]]['wtxnyms'])
  expect_true(trs[[2]]['wspn'])
  expect_true(trs[[2]]['wtxnyms'])
  expect_true(trs[[3]]['ntips'] == 100)
  expect_false(trs[[3]]['wspn'])
  expect_false(trs[[3]]['wtxnyms'])
})
if(file.exists('test.tre')) {
  file.remove('test.tre')
}
if(file.exists('test.trmn')) {
  file.remove('test.trmn')
}