# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'TreeMan Class\'')
test_that('.updateSlots() works', {
  tree <- randTree(100)
  tree@nodelist[['t1']][['span']] <- NULL
  tree <- treeman:::.updateSlots(tree)
  expect_false(tree@wspn)
})
test_that('.localUpdateKids() works', {
  tree <- randTree(10)
  kids <- getNodesKids(tree, tree['nodes'])
  ndlst <- tree['nodelist']
  tid <- tree['tips'][1]
  for(i in 1:length(ndlst)) {
    bool <- ndlst[[i]][['kids']] != tid
    ndlst[[i]][['kids']] <- ndlst[[i]][['kids']][bool]
  }
  ndlst <- treeman:::.localUpdateKids(ndlst, tid=tid,
                                      rid=tree['root'])
  tree@nodelist <- ndlst
  new_kids <- getNodesKids(tree, tree['nodes'])
  res <- rep(NA, length(kids))
  for(i in 1:length(res)) {
    res[i] <- all(kids[[i]] %in% new_kids[[i]]) &
      length(kids[[i]]) == length(new_kids[[i]])
  }
  expect_true(all(res))
})
test_that('.globalUpdateKids() works', {
  tree <- randTree(10)
  kids <- getNodesKids(tree, tree['nodes'])
  ndlst <- tree['nodelist']
  for(i in 1:length(ndlst)) {
    ndlst[[i]][['kids']] <- NULL
  }
  ndlst <- treeman:::.globalUpdateKids(ndlst)
  tree@nodelist <- ndlst
  new_kids <- getNodesKids(tree, tree['nodes'])
  res <- rep(NA, length(kids))
  for(i in 1:length(res)) {
    res[i] <- all(kids[[i]] %in% new_kids[[i]]) &
      length(kids[[i]]) == length(new_kids[[i]])
  }
  expect_true(all(res))
})
test_that('.localUpdateTip() works', {
  tree <- randTree(10)
  ndlst <- tree['nodelist']
  tid <- sample(tree['tips'], 1)
  ndlst <- treeman:::.localDwndateTip(ndlst, tid=tid, rid=tree['root'])
  new_ndlst <- treeman:::.localUpdateTip(ndlst, tid=tid, rid=tree['root'])
  tree <- new('TreeMan', nodelist=new_ndlst, root=tree['root'])
  tree <- treeman:::.updateSlots(tree)
  expect_that(tree['ntips'], equals(10))
})
test_that('.localDwndateTip() works', {
  tree <- randTree(10)
  ndlst <- tree['nodelist']
  tid <- sample(tree['tips'], 1)
  new_ndlst <- treeman:::.localDwndateTip(ndlst, tid=tid, rid=tree['root'])
  pd_res <- kid_res <- rep(NA, length(ndlst))
  for(i in 1:length(ndlst)) {
    kid_res[i] <- !(tid %in% new_ndlst[[i]][['kids']])
    pd_res[i] <- ndlst[[i]][['pd']] > new_ndlst[[i]][['pd']]
  }
  expect_true(all(kid_res))
  expect_that(sum(pd_res), is_more_than(0))
})
test_that('.localDwndateNode() works', {
  tree <- randTree(10)
  ndlst <- tree['nodelist']
  nid <- sample(tree['nodes'][tree['nodes'] != tree['root']], 1)
  new_ndlst <- treeman:::.localDwndateNode(ndlst, nid=nid, rid=tree['root'])
  pd_res <- prdst_res <- rep(NA, length(ndlst))
  for(i in 1:length(ndlst)) {
    pd_res[i] <- ndlst[[i]][['pd']] > new_ndlst[[i]][['pd']]
  }
  expect_that(sum(pd_res), is_more_than(0))
})
test_that('.localUpdateNode() works', {
  tree <- randTree(10)
  ndlst <- tree['nodelist']
  nid <- sample(tree['nodes'][tree['nodes'] != tree['root']], 1)
  ndlst <- treeman:::.localDwndateNode(ndlst, nid=nid, rid=tree['root'])
  new_ndlst <- treeman:::.localUpdateNode(ndlst, nid=nid, rid=tree['root'])
  tree <- new('TreeMan', nodelist=new_ndlst, root=tree['root'])
  tree <- treeman:::.updateSlots(tree)
  expect_that(tree['ntips'], equals(10))
})
test_that('.globalUpdateAll() works', {
  tree <- randTree(10)
  ndlst <- tree['nodelist']
  for(i in 1:length(ndlst)) {
    ndlst[[i]][['prdst']] <- ndlst[[i]][['pd']] <-
      ndlst[[i]][['kids']] <- NULL
  }
  new_ndlst <- treeman:::.globalUpdateAll(ndlst=ndlst)
  tree <- new('TreeMan', nodelist=new_ndlst, root=tree['root'])
  tree <- treeman:::.updateSlots(tree)
  expect_that(tree['ntips'], equals(10))
})