# LIB
library(treeman)

# Load BIG trees and get ages
data(mammals)
data(birds)
mml_age <- getAge(mammals)
brd_age <- getAge(birds)

# Rename internal nodes
mammals['nds'][1:10]
birds['nds'][1:10]  # possible node IDs are not unique
mml_nds <- paste0('mml_', mammals['nds'])
mml_nds[1:10]
avs_nds <- paste0('avs_', birds['nds'])
mammals <- setNdsID(mammals, mammals['nds'], mml_nds)
birds <- setNdsID(birds, birds['nds'], avs_nds)

# Generate incipient tree
tree <- twoer()

# Add to mammals and birds to tree
tree <- setNdsID(tree, c('t1', 't2'), c('Mammalia', 'Aves'))
# set spans to make tree ultrametric
tree <- setNdSpn(tree, 'Aves', abs(brd_age - mml_age) + 1)
# list of form: id=txnym
txnyms <- list('Aves'='Aves', 'Mammalia'= 'Mammalia',
               'root'='Amniota')
tree <- setTxnyms(tree, txnyms)
tree <- addClade(tree=tree, id='Mammalia', clade=mammals)
tree <- addClade(tree=tree, id='Aves', clade=birds)

# Complete
tree[['root']]
summary(tree)
# Note, plotting takes time.... have to convert to APE
plot(as(tree, 'phylo'), show.tip.label=FALSE)
