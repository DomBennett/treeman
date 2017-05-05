# GENERATE DATA TREES
library(treeman)

# MAMMALS
mammals <- readTree(file='other/mammals.tre')
txnyms <- searchTxnyms(mammals, cache=TRUE, infer=TRUE,
                       clean=TRUE, parent='Mammalia')
mammals <- setTxnyms(mammals, txnyms)
summary(mammals)
save(mammals, file="data/mammals.rda", compress="xz")

# BIRDS
birds <- readTree(file='other/birds.tre')
txnyms <- searchTxnyms(birds, cache=TRUE, infer=TRUE,
                       clean=TRUE, parent='Aves')
birds <- setTxnyms(birds, txnyms)
summary(birds)
save(birds, file="data/birds.rda", compress="xz")

# PLANTS
plants <- ape::read.tree(file='other/plants.tre')
plants$root.edge <- NULL
plants$node.label <- NULL  # treeman cannot handle special characters in nodel label
plants <- as(plants, 'TreeMan')
txnyms <- searchTxnyms(plants, cache=TRUE, infer=TRUE,
                       clean=TRUE, parent='Viridiplantae')
plants <- setTxnyms(plants, txnyms)
summary(plants)
save(plants, file="data/plants.rda", compress="xz")

