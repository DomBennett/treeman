# There is an error when reading the large birds tree
# let's read each tree in a loop and see if we can diagnose the problem
library(treeman)
file <- '0_data/trees/birds.tre'
file <- "~/Desktop/birds.tre"
cat('Reading in trstrs...\n', sep="")
trstrs <- scan(file, what="raw", quiet=TRUE)
cat("Found [", length(trstrs), "] strs.\n", sep="")
trees <- vector("list", length=length(trstrs))
cat("Looping ....\n")
if(!file.exists("bird_trees")) {
  dir.create("bird_trees")
}
for(i in 1:length(trstrs)) {
  cat(i, ", ", sep="")
  tree <- treeman:::.readTree(trstrs[i], update=FALSE)
  save(tree, file=file.path("bird_trees", paste0("tree_", i, ".RData")))
  trees[[i]] <- tree
}
cat("\nConverting to TreeMen ....\n")
trees <- as(trees, 'TreeMen')
cat("\nSaving ....\n")
save(trees, file="birds_trees.RData")

eds <- calcFrPrp(tree, tree['tips'])


