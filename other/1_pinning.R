# Pinning mammal species to mammalian supertree using taxonomy

# LIBS
library(treeman)

# DATA
data(mammals)  # example mammal tree is 'taxonomically informed', all nodes have taxonyms
# rslvd names (i.e. with lineages) not in mammals tree generated with taxaResolve()
# load pre-generated dataset from github, requires internet
load(url("https://github.com/DomBennett/treeman/raw/master/other/1_pinning.RData"))

# PARAMETERS
n <- 1000  # number of missing mammal species to pin

# CLEAN DATA
rnds <- sample(1:nrow(rslvd_mammals), n)
rslvd_mammals <- rslvd_mammals[rnds, ]
lngs <- lapply(rslvd_mammals$lineage, function(lineage, ...) strsplit(lineage, '\\|')[[1]])
# remove everything before mammalia
pull <- sapply(lngs, function(lineage) any(grepl('mammalia', lineage, ignore.case=TRUE)))
lngs <- lngs[pull]
lngs <- lapply(lngs, function(lineage) {
  lineage[-1*(1:(which(grepl('mammalia', lineage, ignore.case=TRUE)) - 1))]
  })
# always replace spaces with _
lngs <- lapply(lngs, function(x) gsub("\\s+", "_", x))
tids <- gsub("\\s+", "_", rslvd_mammals$search_name)
tids <- tids[pull]
cat('[', sum(pull)*100/n, '%] random names are pinnable\n', sep='')

# PIN
ends <- rep(0, length(tids))  # all tips end in the present
tree_age <- getAge(mammals)
pinned_tree <- pinTips(tree=mammals, lngs=lngs, tids=tids,
                       end_ages=ends, tree_age=tree_age)
pinned_tree <- updateSlts(pinned_tree)
p_added <- sum(tids %in% pinned_tree['tips'])*100/n
cat('[', p_added, '%] of n pinned to mammals\n', sep='')

# VIZ (using ape)
txnym <- function(n) {
  # return taxonyms for node labels, combining any multiple entries with _
  paste0(n[['txnym']], collapse='_')
}
writeTree(pinned_tree, file='temp.tre', ndLabels = txnym)
tree_phylo <- ape::read.tree('temp.tre')
plot(tree_phylo, show.tip.label=FALSE, edge.width=0.5, type='fan', no.margin=FALSE,
     edge.color='lightsteelblue3')
tids <- tids[tids %in% tree_phylo$tip.label]
# only plot 10 names and only use last two elements of a name -- easier to see
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