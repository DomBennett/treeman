# TESTING FOR SIGNIFICANT PHYLOGENETIC TURNOVER

# LIBS
library(treeman)
# requires two functions from MoreTreeTools
source('other/mtt_community_tools.R')

# PARAMETERS
ntips <- 500
tree <- ape::rtree(ntips)
dpsi <- 5  # power difference

# GENERATE DATA
psi_1 <- dpsi
psi_2 <- -dpsi
focal <- round(ntips*0.5)
mean.incid <- ntips*0.05
c1 <- genCommData(tree=tree, psi=psi_1,
                  mean.incid=mean.incid,
                  mean.abun=mean.incid,
                  nsites=1, focal=focal)
c2 <- genCommData(tree=tree, psi=psi_2,
                  mean.incid=mean.incid,
                  mean.abun=mean.incid,
                  nsites=1, focal=focal)

# VIZ
# construct community matrix for community plot
cmatrix <- rbind(c1, c2)
cmatrix[cmatrix > 0] <- 1
commplot(cmatrix, tree, groups=c(1,2), no.margin=FALSE)
mtext(text=paste0('psi = ', dpsi))

# PERMUTATION TEST WITH TREEMAN
tree_tm <- as(tree, 'TreeMan')
c1_ids <- colnames(c1)[c1[1, ] > 0]
c2_ids <- colnames(c2)[c2[1, ] > 0]
obs_ovrlp <- calcOvrlp(tree_tm, c1_ids, c2_ids)  # determine the proportion of shared branch length
iterations <- 99
null <- rep(NA, iterations)
for(i in 1:iterations) {
  cat('.... [', i, ']\n', sep='')
  null_tips <- sample(tree_tm['tips'], length(c1_ids))  # generate null distributions
  null[i] <- calcOvrlp(tree_tm, c1_ids, null_tips)
}
p_value <- sum(obs_ovrlp >= null)/iterations
hist(null, main="", xlab="", ylab="")
abline(v=obs_ovrlp, col="red")
mtext(paste0("P-value: ", signif(p_value, 3)))
cat("P-value: ", signif(p_value, 3), "\n", sep="")

# Plot both together
par(mfrow=c(1,2))
commplot(cmatrix, tree, groups=c(1,2), no.margin=FALSE)
mtext(text=paste0('psi = ', dpsi))
hist(null, main="", xlab="", ylab="")
abline(v=obs_ovrlp, col="red")
mtext(paste0("P-value: ", signif(p_value, 3)))
