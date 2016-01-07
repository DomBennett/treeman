# TODO: addTips, removeTip, mergeTree, collapseNode
# TODO: add doc for adding and removing tips
addTip <- function(tree, id, sister, start, end,
                   parent_id=paste0("p_", id)) {
  updatePrenodes <- function(node) {
    node$children <- c(node$children, tip$id)
    node$pd <- node$pd + tip$span
    node
  }
  tip <- list('id'=id)
  node <- list('id'=parent_id)
  tip$span <- start - end
  age <- getNodeAge(tree, sister)
  new_sister <- sister <- tree@nodelist[[sister]]
  new_parent <- tree@nodelist[[sister$prenode]]
  new_parent$postnode <- new_parent$postnode[!new_parent$postnode %in% sister$id]
  new_parent$postnode <- c(new_parent$postnode, node$id)
  new_sister$span <- start - age
  new_sister$prenode <- node$id
  node$span <- sister$span - new_sister$span
  node$pd <- new_sister$span + tip$span
  node$predist <- sister$predist - new_sister$span
  node$prenode <- sister$prenode
  node$postnode <- node$children <- c(tip$id, sister$id)
  tip$pd <- 0
  tip$predist <- sister$predist + tip$span
  tip$prenode <- node$id
  tree@nodelist[[tip$id]] <- tip
  tree@nodelist[[node$id]] <- node
  tree@nodelist[[new_sister$id]] <- new_sister
  tree@nodelist[[new_parent$id]] <- new_parent
  prenodes <- getNodePrenodes(tree, node$id)
  tree@nodelist[prenodes] <- lapply(tree@nodelist[prenodes],
                                    updatePrenodes)
  .update(tree)
}
#TODO: add doc on pinning tips
pinTip <- function(tree, tip_id, lineage, end) {
  lineage <- lineage[lineage != tree@root]
  rngs <- getEdgesAge(tree, lineage)
  print("----")
  print(end)
  print(rngs)
  print("----")
  rngs <- rngs[rngs[ ,'max'] > end, ]
  rngs[rngs[ ,'min'] <= end, "min"] <- end
  prbs <- rngs$max - rngs$min
  nd <- as.vector(sample(rngs$edge, prob=prbs, size=1))
  i <- which(rngs$edge == nd)
  start <- runif(min=rngs$min[i], max=rngs$max[i], n=1)
  tree <- addTip(tree, id=tip_id, sister=nd, start=start, end=end)
  tree
}

pinTips <- function(tree, tip_ids, lineages, ends) {
  .pin <- function(i) {
    tree <- pinTip(tree, tip_ids[i], lineages[[i]], ends[i])
    tree <<- tree
  }
  sapply(1:length(tip_ids), .pin)
  tree
}
