# MEASURE RUN DIFFERENCES BETWEEN WNDMTRX=T/F
data(mammals)
summary(mammals)
m1 <- mammals
m2 <- addNdmtrx(mammals)
timings <- data.frame(ndlst=rep(NA, 5),
                      ndmtrx=rep(NA, 5))
rownames(timings) <- c('prids', 'ptids', 'prdsts',
                       'ages', 'frprp')
tree <- m1
for(i in 1:2) {
  timings['prids', i] <- system.time(prids <- getNdsPrids(tree, tree['all']))[[1]]
  timings['ptids', i] <- system.time(ptids <- getNdsPtids(tree, tree['all']))[[1]]
  timings['prdsts', i] <- system.time(prdsts <- getNdsPrdst(tree, tree['all']))[[1]]
  timings['ages', i] <- system.time(ages <- getNdsAge(tree, tree['all'], getAge(tree)))[[1]]
  timings['frprp', i] <- system.time(frprps <- calcFrPrp(tree, tree['tips']))[[1]]
  tree <- m2
}
mean(timings[['ndlst']]/timings[['ndmtrx']]) # ndmtrx is 2-3 times faster
