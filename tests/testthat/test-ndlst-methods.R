# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'ndlst-methods\'')
# slot
test_that('.getSltPrids() works', {
  tree <- randTree(10)
  prids <- treeman:::.getSltPrids(tree@ndlst, FALSE)
  i <- sample(1:tree@nall, 1)
  prid <- tree@ndlst[[i]][['prid']]
  expect_equal(prids[i], prid)
})
test_that('.getSltSpns() works', {
  tree <- randTree(10)
  spns <- treeman:::.getSltSpns(tree@ndlst, FALSE)
  i <- sample(1:tree@nall, 1)
  spn <- tree@ndlst[[i]][['spn']]
  expect_equal(spns[i], spn)
})
test_that('.getSltTids() works', {
  tree <- randTree(10)
  tids <- treeman:::.getSltTids(tree@ndlst, FALSE)
  expect_true(all(tree['tips'] %in% tids))
  expect_true(all(tids %in% tree['tips']))
})
# single node
test_that('.getNdSstrFrmLst() works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  # test
  rnd_tip <- sample(tree@tips, 1)
  res <- treeman:::.getNdSstrFrmLst(tree@ndlst, rnd_tip)
  expect_true(res %in% ids)
})
test_that('.getNdPridsFrmLst() works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  prids <- treeman:::.getSltPrids(tree@ndlst, FALSE)
  # tests
  rnd_tip <- sample(tree@tips, 1)
  res <- treeman:::.getNdPridsFrmLst(ndlst, prids=prids,
                                     id=rnd_tip)
  expect_true(ndlst[[rnd_tip]]['prid'] %in% res)
  expect_true(tree@root %in% res)
  rnd_nd <- sample(tree@nds, 1)
  res <- treeman:::.getNdPridsFrmLst(ndlst, prids=prids,
                                     id=rnd_nd)
  expect_true(ndlst[[rnd_nd]]['prid'] %in% res)
  expect_true(tree@root %in% res)
  res <- treeman:::.getNdPridsFrmLst(ndlst, prids=prids,
                                     tree@root)
  expect_true(tree@root == res)
})
test_that('.getNdPrdstFrmLst([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  prids <- treeman:::.getSltPrids(tree@ndlst, FALSE)
  # tests
  res <- treeman:::.getNdPrdstsFrmLst(ndlst, prids=prids,
                                     id=tree@root)
  expect_that(res, equals(0))
  rnd_nd <- sample(tree@nds[tree@nds != tree@root], 1)
  res <- treeman:::.getNdPrdstsFrmLst(ndlst, prids=prids,
                                     id=rnd_nd)
  expect_that(res, is_less_than(tree@age))
})
test_that('.getNdPtidsFrmLst([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  prids <- treeman:::.getSltPrids(tree@ndlst, FALSE)
  # tests
  rnd_tip <- sample(tree@tips, 1)
  res <- treeman:::.getNdPtidsFrmLst(ndlst, prids=prids,
                                     id=rnd_tip)
  expect_that(length(res), equals(0))
  intrnls <- tree@nds[tree@nds != tree@root]
  res <- treeman:::.getNdPtidsFrmLst(ndlst, prids=prids,
                                     id=tree@root)
  expect_true(all(intrnls %in% res))
  rnd_nd <- sample(intrnls, 1)
  res <- treeman:::.getNdPtidsFrmLst(ndlst, prids=prids,
                                     id=rnd_nd)
  expect_false(all(intrnls %in% res))
})
test_that('.getNdKidsFrmLst([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  prids <- treeman:::.getSltPrids(tree@ndlst, FALSE)
  tids <- treeman:::.getSltTids(tree@ndlst, FALSE)
  # tests
  res <- treeman:::.getNdKidsFrmLst(ndlst, prids=prids,
                                    tids=tids, id=tree@root)
  expect_true(all(tree@tips %in% res))
  rnd_nd <- sample(tree@nds[tree@nds != tree@root], 1)
  res <- treeman:::.getNdKidsFrmLst(ndlst, prids=prids,
                                    tids=tids, id=rnd_nd)
  expect_false(all(tree@tips %in% res))
  rnd_tip <- sample(tree@tips, 1)
  res <- treeman:::.getNdKidsFrmLst(ndlst, prids=prids,
                                    tids=tids, id=rnd_tip)
  expect_that(length(res), equals(0))
})
test_that('.getNdPDFrmLst([basic]) works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  ids <- names(ndlst)
  prids <- treeman:::.getSltPrids(tree@ndlst, FALSE)
  # tests
  res <- treeman:::.getNdPDFrmLst(ndlst, prids=prids,
                                  id=tree@root)
  expect_that(res, equals(tree@pd))
  rnd_nd <- sample(tree@nds[tree@nds != tree@root], 1)
  res <- treeman:::.getNdPDFrmLst(ndlst, prids=prids,
                                  id=rnd_nd)
  expect_that(res, is_less_than(tree@pd))
})
# multiple nodes
test_that('.getNdsPtidsFrmLst() works', {
  #TODO
})
test_that('.getNdsPridsFrmLst() works', {
  #TODO
})
test_that('.getNdsPDFrmLst() works', {
  #TODO
})
test_that('.getNdsKidsFrmLst() works', {
  #TODO
})
test_that('.getNdsPrdstsFrmLst() works', {
  #TODO
})
# tree
test_that('.getTreeAge() works', {
  # init
  tree <- randTree(100)
  ndlst <- tree@ndlst
  # test
  res <- treeman:::.getTreeAgeFrmLst(tree@ndlst, FALSE)
  expect_lt(res, tree['pd'])
  expect_equal(res, tree['age'])
})