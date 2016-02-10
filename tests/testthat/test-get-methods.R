# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'get-methods\'')
test_that('getOutgroup() works', {
  tree <- randTree(10)
  rnd_nd <- sample(tree['nodes'][tree['nodes'] != tree['root']], 1)
  ingrp <- getNodeKids(tree, rnd_nd)
  outgrp <- sample(tree['tips'][!tree['tips'] %in% ingrp], 1)
  res <- getOutgroup(tree, ids=c(ingrp, outgrp))
  expect_that(res, equals(outgrp))
})
test_that('get_Slot() works', {
  tree <- randTree(10)
  node_spans <- getNodesSlot(tree, name="span",
                             ids=tree['all'])
  expect_that(sum(node_spans), equals(tree['pd']))
})
test_that('get_Kids() works', {
  tree <- randTree(10)
  kids <- getNodesKids(tree, tree['nodes'])
  expect_true(all(kids$n1 %in% paste0("t", 1:10)))
})
test_that('get_Age() works', {
  tree <- randTree(10)
  root_age <- tree['age']
  nd_ages <- getNodesAge(tree, tree['nodes'])
  expect_true(all(nd_ages <= root_age))
})
test_that('getSpanAge() works', {
  tree <- randTree(10)
  tip_age <- getSpanAge(tree, sample(tree['tips'], 1))
  expect_that(tip_age['start'], is_more_than(tip_age['end']))
})
test_that('getSpansAge() works', {
  tree <- randTree(10)
  tip_ages <- getSpansAge(tree, tree['tips'])
  res <- all(tip_ages[ ,'start'] > tip_ages[ ,'end'])
  expect_true(res)
})
test_that("getParent() works", {
  tree <- readTree(text="(((A,B),(C,D)),(E,F));")
  prnt <- getParent(tree, ids=c("A", "C"))
  expect_true(prnt == "n2")
})
test_that("getPath() works", {
  tree <- randTree(10)
  pth <- getPath(tree, from="t1", to="t10")
  prnt <- getParent(tree, ids=c('t1', "t10"))
  expect_true(prnt %in% pth)
  expect_that(pth[1], equals('t1'))
  expect_that(pth[length(pth)], equals('t10'))
})
test_that("get_Prid() works", {
  tree <- randTree(10)
  prid <- getNodePrid(tree, id='n1')
  expect_that(prid, is_null())
  prids <- getNodesPrid(tree, tree['nodes'])
  lst_nds <- unlist(lapply(prids, function(n) n[length(n)]))
  expect_true(all(lst_nds == "n1"))
})
test_that("get_Ptid() works", {
  tree <- randTree(10)
  pstids <- getNodesPtid(tree, tree['nodes'])
  n1_ptids <- tree['all'][tree['all'] != 'n1']
  expect_true(all(n1_ptids %in% pstids[['n1']]))
  expect_that(pstids[['t1']], is_null())
})
test_that("get_Lineage() works", {
  # TODO: redo with setNodes
  tree <- randTree(10)
  for(i in 1:length(tree@nodelist)) {
    tree@nodelist[[i]]$taxonym <- paste0("l", sample(1:1000, 1))
  }
  lngs <- getNodesLineage(tree, tree['tips'])
  rnd1 <- sample(1:length(lngs), 1)
  rnd2 <- sample(1:length(lngs), 1)
  expect_that(sum(lngs[[rnd1]] %in% lngs[[rnd2]]),
              is_more_than(0))
})
test_that("getSubtree() works", {
  tree <- randTree(10)
  subtree <- getSubtree(tree, 'n2')
  expect_that(tree['ntips'], is_more_than(subtree['ntips']))
  expect_that(tree['nnodes'], is_more_than(subtree['nnodes']))
  expect_that(tree['pd'], is_more_than(subtree['pd']))
  expect_that(tree['age'], is_more_than(subtree['age']))
})