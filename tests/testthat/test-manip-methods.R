# LIBS
library(treeman)
library(testthat)

# TEST FUNCTIONS
randomTips <- function(n, tree) {
  # generate random tips for pinning
  lineages <- ends <- tip_ids <- rep(NA, n)
  nodes <- names(tree@nodelist)
  nodes <- nodes[nodes != 'n1']
  for (i in 1:n) {
    random_node <- sample(nodes, 1)
    lineages[i] <- list(c(getNodePrenodes(tree, random_node), random_node))
                        ends[i] <- runif(max=tree@age, min=0, n=1)
                        tip_ids[i] <- paste0('new_', i)
  }
  lineages <<- lineages
  ends <<- ends
  tip_ids <<- tip_ids
}

# RUNNING
context('Testing \'manip-methods\'')
test_that('addTip() works', {
  # random tree + basic stats
  tree <- randTree(10)
  pd_before <- tree@pd
  age_before <- tree@age
  ntips_before <- nTips(tree)
  # add random tip
  sister <- sample(tree@tips, 1)
  sister_age <- getNodeAge(tree, sister)
  parent_age <- getNodeAge(tree, tree@nodelist[[sister]]$prenode)
  start <- runif(min=sister_age, max=parent_age, n=1)
  end <- runif(min=0, max=start, n=1)
  tree <- addTip(tree, id='new_tip', sister=sister, start=start, end=end)
  #viz(tree)
  # test if successful
  expect_that(tree@plytms, is_false())
  expect_that(tree@age, equals(age_before))
  expect_that(nTips(tree), equals(ntips_before + 1))
  expect_that(tree@pd, equals(pd_before + (start-end)))
})

test_that('pinTip() and pinTips() work', {
  tree <- randTree(5)
  pd_before <- tree@pd
  age_before <- tree@age
  lineages <- ends <- tip_ids <- NULL
  randomTips(2, tree)
  tree <- pinTips(tree, tip_ids, lineages, ends)
  test_that(nTips(tree), equals(20))
  test_that(pd_before, is_less_than(tree@pd))
  test_that(age_before, equals(tree@age))
})

viz(tree)
tree_1 <- pinTip(tree, tip_id=tip_ids[1], lineage=lineages[[1]],
                 end=ends[1])
