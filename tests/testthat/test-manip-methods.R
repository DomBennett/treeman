# LIBS
library(treeman)
library(testthat)

# TEST FUNCTIONS
randomLineage <- function(n, tree) {
  # add random monophyletic taxonyms to tree
  addLname <- function(node, tree) {
    tree@nodelist[[node]][['txnym']] <- lname
    pnodes <- tree@nodelist[[node]][['ptid']]
    if(!is.null(pnodes)) {
      for(pnode in pnodes) {
        tree <- addLname(pnode, tree)
      }
    }
    tree
  }
  nodes <- tree@nodes
  nodes <- nodes[nodes != tree@root]
  lname <- paste0('l', 1)
  tree <- addLname(tree@root, tree)
  for(i in 2:(n-1)) {
    lname <- paste0('l', i)
    node <- sample(nodes, 1)
    tree <- addLname(node, tree)
  }
  tree
}
randomTips <- function(n, tree) {
  # generate random tips with lineages for pinning
  lineages <- ends <- tip_ids <- rep(NA, n)
  nodes <- names(tree@nodelist)
  nodes <- nodes[nodes != 'n1']
  for (i in 1:n) {
    random_node <- sample(nodes, 1)
    l <- c(getNodeLineage(tree, random_node),
           paste0('new_l', i))
    lineages[i] <- list(l)
    ends[i] <- runif(max=tree@age, min=0, n=1)
    tip_ids[i] <- paste0('new_', i)
  }
  list("l"=lineages, "e"=ends, "t"=tip_ids)
}

# RUNNING
context('Testing \'manip-methods\'')
test_that('addTip() works', {
  # random tree + basic stats
  tree <- randTree(10)
  pd_before <- tree['pd']
  age_before <- tree['age']
  ntips_before <- tree['ntips']
  # add random tip
  sister <- sample(tree@tips, 1)
  sister_age <- getNodeAge(tree, sister)
  parent_age <- getNodeAge(tree, tree@nodelist[[sister]][['prid']][1])
  start <- runif(min=sister_age, max=parent_age, n=1)
  end <- runif(min=0, max=start, n=1)
  tree <- addTip(tree, id='new_tip', sister=sister, start=start, end=end)
  #viz(tree)
  # test if successful
  expect_that(validObject(tree), is_true())
  expect_that(tree['ply'], is_false())
  expect_that(tree['age'], equals(age_before))
  expect_that(tree['ntips'], equals(ntips_before + 1))
  expect_that(tree['pd'], equals(pd_before + (start-end)))
})

test_that('pinTip() and pinTips() work', {
  n_start <- 25
  n_add <- 5
  tree <- randTree(n_start)
  tree <- randomLineage(n_start/2, tree)
  pd_before <- tree['pd']
  age_before <- tree['age']
  rdata <- randomTips(n_add, tree)
  tree <- pinTips(tree, tids=rdata[["t"]],
                  lngs=rdata[["l"]],
                  ends=rdata[["e"]])
  writeTree(tree, file='test.tre')  # expect no error
  expect_that(validObject(tree), is_true())
  #expect_that(tree['ntips'], equals(n_start+n_add))
  expect_that(pd_before, is_less_than(tree['pd']))
  expect_that(age_before, equals(tree['age']))
})
if(file.exists('test.tre')) {
  file.remove('test.tre')
}