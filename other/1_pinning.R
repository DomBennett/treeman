# Pinning mammal species to mammalian supertree using taxonomy

# LIBS
library(treemantools)  # package with online resolution functions

# DATA
data(mammals)
data(rslvd_mammals)  # rslvd names (i.e. with lineages) not in mammals tree

# PARAMETERS
n <- 100  # number of missing mammal species to pin

# PIN
rnds <- sample(1:nrow(rslvd_mammals), n)
rslvd_mammals <- rslvd_mammals[rnds, ]
lngs <- mlply(rslvd_mammals, function(lineage, ...) strsplit(lineage, '\\|')[[1]])
ids <- gsub("\\s+", "_", rslvd_mammals$search_name)
ends <- rep(0, length(ids))  # all tips end in the present
pinned_tree <- pinTips(tree=mammals, lngs=lngs, tids=ids, ends=ends)
ids %in% pinned_tree['tips']

# VIZ
tree_phylo <- as(mammals, 'phylo')
plot(tree_phylo)
