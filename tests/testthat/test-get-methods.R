# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'get-methods\'')
test_that('get_Children() works', {
  tree <- randTree(10)
  children <- getNodesChildren(tree, nodes(tree))
  expect_true(all(children$n1 %in% paste0("t", 1:10)))
})
test_that('get_Age() works', {
  tree <- randTree(10)
  root_age <- age(tree)
  nd_ages <- getNodesAge(tree, nodes(tree))[,2]
  expect_true(all(nd_ages <= root_age))
})
test_that("getParent() works", {
  tree <- readTree(tree_string="(((A,B),(C,D)),(E,F));")
  prnt <- getParent(tree, nodes=c("A", "C"))
  expect_true(prnt == "n2")
})
test_that("getPath() works", {
  tree <- randTree(10)
  pth <- getPath(tree, from="t1", to="t10")
  prnt <- getParent(tree, nodes=c('t1', "t10"))
  expect_true(prnt %in% pth)
})
test_that("get_Prenodes() works", {
  tree <- randTree(10)
  prnds <- getNodesPrenodes(tree, nodes(tree))
  expect_that(prnds[[1]], is_null())
  lst_nds <- unlist(lapply(prnds[-1], function(n) n[length(n)]))
  expect_true(all(lst_nds == "n1"))
})
test_that("get_Postnodes() works", {
  tree <- randTree(10)
  pstnds <- getNodesPostnodes(tree, nodes(tree))
  expect_true(all(pstnds[['n1']] %in% c(nodes(tree), tips(tree))))
  expect_that(pstnds[['t1']], is_null())
})
test_that("get_Lineage() works", {
  # TODO: redo with setNodes
  tree <- randTree(10)
  for(i in 1:length(tree@nodelist)) {
    tree@nodelist[[i]]$taxonym <- paste0("l", sample(1:1000, 1))
  }
  lngs <- getNodesLineage(tree, tips(tree))
  rnd1 <- sample(1:length(lngs), 1)
  rnd2 <- sample(1:length(lngs), 1)
  expect_that(sum(lngs[[rnd1]] %in% lngs[[rnd2]]),
              is_more_than(0))
})
test_that("getSubtree() works", {
  tree <- randTree(10)
  subtree <- getSubtree(tree, 'n2')
  expect_that(nTips(tree), is_more_than(nTips(subtree)))
  expect_that(nNodes(tree), is_more_than(nNodes(subtree)))
  expect_that(pd(tree), is_more_than(pd(subtree)))
  expect_that(age(tree), is_more_than(age(subtree)))
})