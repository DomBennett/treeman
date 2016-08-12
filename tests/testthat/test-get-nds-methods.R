# LIBS
library(treeman)
library(testthat)

# PARAMETERS
n <- 10  # number of tips in test trees

# FUNCTIONS
getTestTree <- function(n) {
  # test 'update'
  update <- sample(c(TRUE, FALSE), 1)
  randTree(n, update)
}

# RUNNING
context('Testing \'get-methods\'')
test_that('getNdsPD() works', {
  tree <- getTestTree(n)
  tot_pd <- sum(getNdsSlt(tree, 'spn', tree['all']))
  pds <- getNdsPD(tree, tree['all'])
  expect_true(all(pds[tree['tips']] == 0))
  expect_true(all(pds[tree['tips']] <= tot_pd))
})
test_that('getNdsPrdst() works', {
  tree <- getTestTree(n)
  ids <- tree['all']
  tree_age <- getTreeAge(tree)
  prdsts <- getNdsPrdst(tree, ids)
  expect_true(all(prdsts <= tree_age))
})
test_that('getNdsSstr() works', {
  tree <- getTestTree(n)
  tips <- tree['tips']
  fwrd <- getNdsSstr(tree, tips)
  rvrse <- getNdsSstr(tree, fwrd)
  expect_that(rvrse, equals(tips))
})
test_that('getNdsSlt() works', {
  tree <- getTestTree(n)
  tot_pd <- sum(getNdsSlt(tree, 'spn', tree['all']))
  nd_spns <- getNdsSlt(tree, slt_nm="spn",
                       ids=tree['all'])
  expect_that(sum(nd_spns), equals(tot_pd))
})
test_that('getNdsKids() works', {
  tree <- getTestTree(n)
  kids <- getNdsKids(tree, tree['nds'])
  expect_true(all(kids$n1 %in% paste0("t", 1:10)))
})
test_that('getNdsAge() works', {
  tree <- getTestTree(n)
  root_age <- tree['age']
  nd_ages <- getNdsAge(tree, tree['nds'], tree_age=tree['age'])
  expect_true(all(nd_ages <= root_age))
})
test_that('getSpnsAge() works', {
  tree <- getTestTree(n)
  tree_age <- getTreeAge(tree)
  tip_ages <- getSpnsAge(tree, tree['tips'], tree_age=tree_age)
  res <- all(tip_ages[ ,'start'] > tip_ages[ ,'end'])
  expect_true(res)
})
test_that("getNdsPrids() works", {
  tree <- getTestTree(n)
  prids <- getNdsPrids(tree, tree['nds'])
  tests <- sapply(prids, function(n) 'n1' %in% n)
  expect_true(all(tests))
  # ordrd
  prids <- getNdsPrids(tree, tree['nds'], ordrd=TRUE)
  lst_ns <- sapply(prids, function(n) n[length(n)])
  expect_true(all(lst_ns == "n1"))
})
test_that("getNdsPtids() works", {
  tree <- getTestTree(n)
  ptids <- getNdsPtids(tree, tree['nds'])
  n1_ptids <- tree['all'][tree['all'] != 'n1']
  expect_true(all(n1_ptids %in% ptids[['n1']]))
  expect_that(ptids[['t1']], is_null())
})
# test_that("get_Lng() works", {
#   # TODO: redo with setNds
#   tree <- getTestTree(n)
#   for(i in 1:length(tree@ndlst)) {
#     tree@ndlst[[i]]$txnym <- paste0("l", sample(1:1000, 1))
#   }
#   lngs <- getNdsLng(tree, tree['tips'])
#   rnd1 <- sample(1:length(lngs), 1)
#   rnd2 <- sample(1:length(lngs), 1)
#   expect_that(sum(lngs[[rnd1]] %in% lngs[[rnd2]]),
#               is_more_than(0))
# })
# test_that('getTxnyms() works', {
#   data('mammals')
#   nid <- sample(mammals['nds'], 1)
#   txnym <- mammals[[nid]]['txnym']
#   res <- getTxnyms(mammals, txnym)
#   expect_true(nid %in% res[[1]])
# })