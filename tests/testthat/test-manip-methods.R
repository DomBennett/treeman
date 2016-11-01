# LIBS
library(treeman)
library(testthat)

# TEST FUNCTIONS
randomLineage <- function(n, tree) {
  # add random monophyletic taxonyms to tree
  addLname <- function(nd, tree) {
    tree@ndlst[[nd]][['txnym']] <- lname
    pnds <- getNdPtids(tree, nd)
    for(pnd in pnds) {
      tree <- addLname(pnd, tree)
    }
    tree
  }
  nds <- tree@nds
  nds <- nds[nds != tree@root]
  lname <- paste0('l', 1)
  tree <- addLname(tree@root, tree)
  for(i in 2:(n-1)) {
    lname <- paste0('l', i)
    nd <- sample(nds, 1)
    tree <- addLname(nd, tree)
  }
  getNdsLng(tree, tree['tips'])
  tree
}
randomTips <- function(n, tree) {
  # generate random tips with lineages for pinning
  lngs <- ends <- tip_ids <- rep(NA, n)
  nds <- names(tree@ndlst)
  nds <- nds[nds != 'n1']
  for (i in 1:n) {
    random_nd <- sample(nds, 1)
    l <- c(getNdLng(tree, random_nd),
           paste0('new_l', i))
    lngs[i] <- list(l)
    ends[i] <- runif(max=tree@age, min=0, n=1)
    tip_ids[i] <- paste0('new_', i)
  }
  list("l"=lngs, "e"=ends, "t"=tip_ids)
}

# RUNNING
context('Testing \'manip-methods\'')
test_that('addTip() works', {
  test_tree_size <- 10
  # TEST 1 check for tree without spns
  tree <- randTree(test_tree_size)
  tree <- setNdsSpn(tree, tree['all'], 0)
  tree <- updateTree(tree)
  tid <- paste0('t', tree['ntips'] + 1)
  sid <- sample(tree['tips'], 1)
  new_tree <- addTip(tree, tid=tid, sid=sid)
  new_tree <- updateTree(new_tree)
  # test if successful
  expect_that(length(getNdKids(new_tree, tree['root'])),
              equals(new_tree['ntips']))
  # TEST 2 check for multiple tips with spns
  tree <- randTree(test_tree_size)
  pd_before <- tree['pd']
  age_before <- tree['age']
  tid <- paste0('t', tree['ntips'] + 1)
  sid <- sample(tree['tips'], 1, replace=FALSE)
  # create suitable age range
  sid_spn <- getNdSlt(tree, 'spn', sid)
  min_strt_age <- getNdAge(tree, sid, tree['age'])
  max_strt_age <- min_strt_age + sid_spn
  strt_age <- runif(1, max=max_strt_age, min=min_strt_age)
  end_age <- strt_age - runif(1, max=strt_age, min=0)
  additional_pd <- sum(strt_age - end_age)
  new_tree <- addTip(tree, tid=tid, sid=sid, strt_age=strt_age,
                      end_age=end_age)
  new_tree <- updateTree(new_tree)
  # phylo <- as(new_tree, 'phylo')
  # plot(phylo)
  spns <- getNdsSlt(tree, 'spn', tree['all'])
  expect_true(all(spns >= 0))
  expect_that(new_tree['age'], equals(age_before))
  expect_that(new_tree['ntips'], equals(test_tree_size + 1))
  expect_that(new_tree['pd'], equals(pd_before + additional_pd))
})
test_that('rmTips() work', {
  n <- 100
  tree <- randTree(n)
  # test dropping 1 tip
  pd_before <- tree['pd']
  tid <- sample(tree['tips'], 1)
  tid_spn <- getNdSlt(tree, id=tid, slt_nm='spn')
  tree <- rmTips(tree, tid)
  expect_true(checkTreeMan(tree))
  expect_that(tree['ntips'], equals(n - 1))
  expect_that(pd_before-tid_spn, equals(tree['pd']))
  # test multiple tips
  n <- tree['ntips']
  pd_before <- tree['pd']
  tids <- sample(tree['tips'], 10)
  tree <- rmTips(tree, tids)
  expect_true(checkTreeMan(tree))
  expect_that(tree['ntips'], equals(n - 10))
  expect_that(pd_before, is_more_than(tree['pd']))
})
test_that('pinTips() work', {
  n_start <- 10
  n_add <- 20
  tree <- randTree(n_start)
  tree <- randomLineage(n_start/2, tree)
  tree <- updateTree(tree)
  pd_before <- tree['pd']
  age_before <- tree['age']
  rdata <- randomTips(n_add, tree)
  tree <- pinTips(tree, tids=rdata[["t"]],
                  lngs=rdata[["l"]],
                  end_ages=rdata[["e"]],
                  tree_age=tree['age'])
  tree <- updateTree(tree)
  # phylo <- as(tree, 'phylo')
  # plot(phylo)
  expect_that(validObject(tree), is_true())
  #expect_that(tree['ntips'], equals(n_start+n_add))  # not necessarily true
  expect_that(pd_before, is_less_than(tree['pd']))
  expect_that(tree[['new_1']]['txnym'], is_a('character'))
  #expect_that(age_before, equals(tree['age']))  # not necessarily true
  writeTree(tree, file='test.tre')  # expect no error
})
if(file.exists('test.tre')) {
  file.remove('test.tre')
}