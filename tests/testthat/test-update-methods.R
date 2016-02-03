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
test_that('.localUpdateChildren() works', {
  tree <- randTree(10)
  children <- getNodesChildren(tree, tree['nodes'])
  ndlst <- tree['nodelist']
  tid <- tree['tips'][1]
  for(i in 1:length(ndlst)) {
    bool <- ndlst[[i]][['children']] != tid
    ndlst[[i]][['children']] <- ndlst[[i]][['children']][bool]
  }
  ndlst <- treeman:::.localUpdateChildren(ndlst, tid=tid,
                                          rid=tree['root'])
  tree@nodelist <- ndlst
  new_children <- getNodesChildren(tree, tree['nodes'])
  res <- rep(NA, length(children))
  for(i in 1:length(res)) {
    res[i] <- all(children[[i]] %in% new_children[[i]]) &
      length(children[[i]]) == length(new_children[[i]])
  }
  expect_true(all(res))
})
test_that('.globalUpdateChildren() works', {
  tree <- randTree(10)
  children <- getNodesChildren(tree, tree['nodes'])
  ndlst <- tree['nodelist']
  for(i in 1:length(ndlst)) {
    ndlst[[i]][['children']] <- NULL
  }
  ndlst <- treeman:::.globalUpdateChildren(ndlst, tids=tree['tips'], rid=tree['root'])
  tree@nodelist <- ndlst
  new_children <- getNodesChildren(tree, tree['nodes'])
  res <- rep(NA, length(children))
  for(i in 1:length(res)) {
    res[i] <- all(children[[i]] %in% new_children[[i]]) &
      length(children[[i]]) == length(new_children[[i]])
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
    kid_res[i] <- !(tid %in% new_ndlst[[i]][['children']])
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
      ndlst[[i]][['children']] <- NULL
  }
  nids <- tree['nodes'][tree['nodes'] != tree['root']]
  new_ndlst <- treeman:::.globalUpdateAll(ndlst=ndlst, tids=tree['tips'],
                                nids=nids, rid=tree['root'])
  tree <- new('TreeMan', nodelist=new_ndlst, root=tree['root'])
  tree <- treeman:::.updateSlots(tree)
  expect_that(tree['ntips'], equals(10))
})