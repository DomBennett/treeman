# Work out memory use for different tree sizes

# LIBS
library(treeman)
library(ggplot2)

# CALC TREE SIZES
tree_n <- round((2:21)^exp(1))  # log scale
tree_sizes <- rep(NA, length(tree_n))
for(i in 1:length(tree_n)) {
  cat(i, " ")
  tree <- randTree(tree_n[i])
  tree_sizes[i] <- as.numeric(object.size(tree))
}

# PLOT
pdata <- data.frame(n=tree_n, mb=tree_sizes/10^6)
p <- ggplot(pdata, aes(x=n, y=mb)) +
  geom_point(colour="cornflowerblue") +
  theme_bw() + ylab("Memory (Mb)") + xlab("N tips")
jpeg(file.path('other', 'tree_memory.jpg'))
print(p)
dev.off()

