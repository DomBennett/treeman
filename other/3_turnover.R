# TESTING FOR SIGNIFICANT PHYLOGENETIC TURNOVER

# LIBS
library(treeman)

# INPUT
# load tree, c1 (community 1) and c2 (community 2) from treeman/other
load(file.path('other', '3_turnover.Rd'))
# calcPhyDv(tree, ids=c1)  # both have comparable PDs
# calcPhyDv(tree, ids=c2)

# RUN
obs_ovrlp <- calcOvrlp(tree, c1, c2)  # determine the proportion of shared branch length
iterations <- 999
null <- rep(NA, iterations)
for(i in 1:iterations) {
  print(i)
  null_tips <- sample(tree['tips'], length(c1))  # generate null distributions
  null[i] <- calcOvrlp(tree, c1, null_tips)
}
p_value <- sum(obs_ovrlp >= null)/iterations

# VIZ
# load tree viz libraries
library(MoreTreeTools)
library(treemantools)
# convert to phylo
tree_phylo <- as(tree, 'phylo')
# construct community matrix for community plot
cmatrix <- matrix(rep(0, tree['ntips']*2), nrow=2)
colnames(cmatrix) <- tree['tips']
cmatrix[1, c1] <- 1
cmatrix[2, c2] <- 1
commplot(cmatrix, tree_phylo, groups=c(1,2))
hist(null, main="", xlab="Shared Phylogenetic Diversity")
abline(v=obs_ovrlp, col="red")
text(x=obs_ovrlp, y=170, labels="Observed", cex=0.75, pos=4, col="red")
cat("P-value: ", signif(p_value, 3), "\n", sep="")