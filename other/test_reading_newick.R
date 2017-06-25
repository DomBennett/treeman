# Testing readTree
# Reading trees downloaded from TreeBase

library(treeman)
trdir <- file.path('other', 'newick_examples')
# tree information
# > 2,000 rooted trees from https://doi.org/10.1017/pab.2016.36
# these trees can be read in and plotted by APE
treeinfo <- read.csv(file.path(trdir, 'treeinfo.csv'))
colnames(treeinfo)
# loop through files
trfls <- list.files(trdir, pattern='.tre')
read_failed <- write_failed <- NULL
read_cntr <- write_cntr <- 0
for(trfl in trfls[1:10]) {
  try(expr={
    tree <- readTree(file=file.path(trdir, trfl))
    writeTree(tree, file='temp.tre')
    },
      silent=TRUE)
  if(exists('tree')) {
    read_cntr <- read_cntr + 1
  } else {
    read_failed <- c(read_failed, trfl)
  }
  suppressWarnings(rm(tree))
  if(file.exists('temp.tre')) {
    write_cntr <- write_cntr + 1
  } else {
    write_failed <- c(write_failed, trfl)
  }
  suppressWarnings(file.remove('temp.tre'))
}
cat('Read success rate: [', signif(read_cntr/length(trfls), digits=4),
    ']\n', sep='')
# >99% success rate, failed for certain bad trees in folder
cat('Write success rate: [', signif(write_cntr/length(trfls), digits=4),
    ']\n', sep='')
