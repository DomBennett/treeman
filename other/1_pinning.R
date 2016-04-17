# Pinning mammal species to mammalian supertree using taxonomy

# LIBS
library(treeman)

# DATA
data(mammals)  # example mammal tree is 'taxonomically informed', all nodes have taxonyms
# rslvd names (i.e. with lineages) not in mammals tree generated with MoreTreeTools::taxaResolve
# load pre-generated dataset from github
load(url("https://github.com/DomBennett/treeman/raw/master/other/1_pinning.RData"))

# PARAMETERS
n <- 100  # number of missing mammal species to pin

# PIN
rnds <- sample(1:nrow(rslvd_mammals), n)
rslvd_mammals <- rslvd_mammals[rnds, ]
lngs <- plyr::mlply(rslvd_mammals, function(lineage, ...) strsplit(lineage, '\\|')[[1]])
tids <- gsub("\\s+", "_", rslvd_mammals$search_name)  # always replace spaces with _
ends <- rep(0, length(tids))  # all tips end in the present
pinned_tree <- pinTips(tree=mammals, lngs=lngs, tids=tids, ends=ends)
p_added <- sum(tids %in% pinned_tree['tips'])*100/n
cat('[', p_added, '%] of n pinned to mammals\n', sep='')

# VIZ
library(MoreTreeTools)  # for conversion to phylo
txnym <- function(n) {
  # return taxonyms for node labels, combining any multiple entries with _
  paste0(n[['txnym']], collapse='_')
}
writeTree(pinned_tree, file='temp.tre', ndLabels = txnym)
tree_phylo <- ape::read.tree('temp.tre')
plot(tree_phylo, show.tip.label=FALSE, edge.width=0.5, type='fan', no.margin=FALSE,
     edge.color='lightsteelblue3')
tids <- tids[tids %in% tree_phylo$tip.label]
# only plot 10 names and only use last two elements of a names -- easier to see
print_tids <- sample(tids, 10, prob=1/sapply(tids, nchar))
print_names <- unlist(lapply(strsplit(print_tids, '_'),
                             function(x) paste0(x[(length(x)-1):length(x)], collapse=' ')))
indx <- match(print_tids, tree_phylo$tip.label)
ape::tiplabels(pch=19, tip=indx, cex=0.9, col='black')
ape::tiplabels(text=print_names,
               tip=indx, cex=0.85, adj=c(-0.05, -0.25),
               bg=NULL, col='black', frame='none', font=3)
# example large clades for guidance
node_labels <- c('Eutheria', 'Metatheria', 'Primates',
                 'Rodentia', 'Chiroptera', 'Carnivora',
                 'Afrotheria', 'Cetartiodactyla')
nds <- match(node_labels, tree_phylo$node.label)+length(tree_phylo$tip.label)
ape::nodelabels(pch=19, node=nds, cex=0.9, col='gray40')
ape::nodelabels(text=node_labels, node=nds, cex=0.85, adj=-0.1,
                bg=NULL, col='gray40', frame='none')
if(file.exists('temp.tre')) {
  file.remove('temp.tre')
}