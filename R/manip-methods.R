# TODO: removeTip, mergeTree, collapseNode

addTip <- function(tree, id, sister, start, end) {
  updatePrenodes <- function(node) {
    node$children <- c(node$children, tip$id)
    node$pd <- node$pd + tip$span
    node
  }
  tip <- list('id'=id)
  node <- list('id'=paste0("p", id))
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