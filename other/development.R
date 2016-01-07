# TODO:
# -- read/write newick

# TODO -- why does this lead to an increase in age sometimes?
test_tree <- randTree(10)
rng <- getEdgeAge(test_tree, 'n8')
start <- runif(max=rng$max, min=rng$min, n=1)
end <- runif(max=start, min=0, n=1)
test_tree2 <- addTip(test_tree, id='new', sister='n8', end=end, start=start)
test_tree2@age == test_tree@age
