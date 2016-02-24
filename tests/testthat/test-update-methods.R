# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'TreeMan Class\'')
test_that('.updateTreeSlots() works', {
  tree <- randTree(100)
  tree@nodelist[['t1']][['span']] <- NULL
  tree <- treeman:::.updateTreeSlots(tree)
  expect_false(tree@wspn)
})
test_that('.updateKids() works', {
  tree <- randTree(10)
  kids <- getNodesKids(tree, tree['nodes'])
  ndlst <- tree['nodelist']
  tid <- sample(tree['tips'], 1)
  for(i in 1:length(ndlst)) {
    bool <- ndlst[[i]][['kids']] != tid
    ndlst[[i]][['kids']] <- ndlst[[i]][['kids']][bool]
  }
  ndlst <- treeman:::.updateKids(ndlst, tid=tid,
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
test_that('.dwndateKids() works', {
  tree <- randTree(10)
  ndlst <- tree['nodelist']
  tid <- sample(tree['tips'], 1)
  ndlst <- treeman:::.dwndateKids(ndlst, tid, tree['root'])
  test_bool <- rep(NA, length(ndlst))
  for(i in 1:length(ndlst)) {
    test_bool[i] <- tid %in% ndlst[[i]][['kids']]
  }
  expect_false(any(test_bool))
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
test_that('.updateTip() works', {
  tree <- randTree(10)
  ndlst <- tree['nodelist']
  tid <- sample(tree['tips'], 1)
  ndlst <- treeman:::.dwndateTip(ndlst, tid=tid, rid=tree['root'])
  new_ndlst <- treeman:::.updateTip(ndlst, tid=tid, rid=tree['root'])
  tree <- new('TreeMan', nodelist=new_ndlst, root=tree['root'])
  tree <- treeman:::.updateTreeSlots(tree)
  expect_that(tree['ntips'], equals(10))
})
test_that('.dwndateTip() works', {
  tree <- randTree(10)
  ndlst <- tree['nodelist']
  tid <- sample(tree['tips'], 1)
  new_ndlst <- treeman:::.dwndateTip(ndlst, tid=tid, rid=tree['root'])
  pd_res <- kid_res <- rep(NA, length(ndlst))
  for(i in 1:length(ndlst)) {
    kid_res[i] <- !(tid %in% new_ndlst[[i]][['kids']])
    pd_res[i] <- ndlst[[i]][['pd']] > new_ndlst[[i]][['pd']]
  }
  expect_true(all(kid_res))
  expect_that(sum(pd_res), is_more_than(0))
})
test_that('.dwndateNode() works', {
  tree <- randTree(10)
  ndlst <- tree['nodelist']
  nid <- sample(tree['nodes'][tree['nodes'] != tree['root']], 1)
  new_ndlst <- treeman:::.dwndateNode(ndlst, nid=nid, rid=tree['root'])
  pd_res <- prdst_res <- rep(NA, length(ndlst))
  for(i in 1:length(ndlst)) {
    pd_res[i] <- ndlst[[i]][['pd']] > new_ndlst[[i]][['pd']]
  }
  expect_that(sum(pd_res), is_more_than(0))
})
test_that('.updateNode() works', {
  tree <- randTree(10)
  age_before <- tree['age']
  pd_before <- tree['pd']
  ntips_before <- tree['ntips']
  ndlst <- tree['nodelist']
  nid <- sample(tree['nodes'][tree['nodes'] != tree['root']], 1)
  ndlst <- treeman:::.dwndateNode(ndlst, nid=nid, rid=tree['root'])
  new_ndlst <- treeman:::.updateNode(ndlst, nid=nid, rid=tree['root'])
  tree <- new('TreeMan', nodelist=new_ndlst, root=tree['root'])
  tree <- treeman:::.updateTreeSlots(tree)
  expect_that(tree['ntips'], equals(ntips_before))
  expect_that(tree['age'], equals(age_before))
  expect_that(tree['pd'], equals(pd_before))
})
test_that('.globalUpdateAll() works', {
  tree <- randTree(10)
  ndlst <- tree['nodelist']
  for(i in 1:length(ndlst)) {
    ndlst[[i]][['prdst']] <- ndlst[[i]][['kids']] <- NULL
    ndlst[[i]][['pd']] <- 0
  }
  new_ndlst <- treeman:::.globalUpdateAll(ndlst=ndlst)
  tree <- new('TreeMan', nodelist=new_ndlst, root=tree['root'])
  tree <- treeman:::.updateTreeSlots(tree)
  expect_that(tree['ntips'], equals(10))
})