# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'get-methods\'')
test_that('get_Slot () works', {
  tree <- randTree(10)
  node_spans <- getNodesSlot(tree, name="span",
                             ids=tree['all'])
  expect_that(sum(node_spans), equals(tree['pd']))
})
test_that('get_Children() works', {
  tree <- randTree(10)
  children <- getNodesChildren(tree, tree['nodes'])
  expect_true(all(children$n1 %in% paste0("t", 1:10)))
})
test_that('get_Age() works', {
  tree <- randTree(10)
  root_age <- tree['age']
  nd_ages <- getNodesAge(tree, tree['nodes'])[,2]
  expect_true(all(nd_ages <= root_age))
})
test_that("getParent() works", {
  tree <- readTree(tree_string="(((A,B),(C,D)),(E,F));")
  prnt <- getParent(tree, ids=c("A", "C"))
  expect_true(prnt == "n2")
})
test_that("getPath() works", {
  tree <- randTree(10)
  pth <- getPath(tree, from="t1", to="t10")
  prnt <- getParent(tree, ids=c('t1', "t10"))
  expect_true(prnt %in% pth)
})
test_that("get_Prid() works", {
  tree <- randTree(10)
  prids <- getNodesPrid(tree, tree['nodes'])
  expect_that(prids[[1]], is_null())
  lst_nds <- unlist(lapply(prids[-1], function(n) n[length(n)]))
  expect_true(all(lst_nds == "n1"))
})
test_that("get_Ptid() works", {
  tree <- randTree(10)
  pstids <- getNodesPtid(tree, tree['nodes'])
  expect_true(all(pstids[['n1']] %in% tree['all']))
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