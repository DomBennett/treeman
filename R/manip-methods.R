# TODO: addTips, removeTip, mergeTree, collapseNode, removeNode
# TODO: add doc for adding and removing tips
addTip <- function(tree, id, sister, start, end,
                   parent_id=paste0("p_", id),
                   tip_taxonym=NULL, parent_taxonym=NULL) {
  updatePre <- function(node) {
    node[['children']] <- c(node[['children']], tip[['id']])
    node[['pd']] <- node[['pd']] + tip[['span']]
    node
  }
  tip <- list('id'=id)
  if(!is.null(tip_taxonym)) {
    tip[['taxonym']] <- tip_taxonym
  }
  node <- list('id'=parent_id)
  if(!is.null(parent_taxonym)) {
    node[['taxonym']] <- parent_taxonym
  }
  tip[['span']] <- start - end
  age <- getNodeAge(tree, sister)
  new_sister <- sister <- tree@nodelist[[sister]]
  new_parent <- tree@nodelist[[sister[['prid']]]]
  new_parent[['ptid']] <- new_parent[['ptid']][!new_parent[['ptid']] %in% sister[['id']]]
  new_parent[['ptid']] <- c(new_parent[['ptid']], node[['id']])
  new_sister[['span']] <- start - age
  new_sister[['prid']] <- node[['id']]
  node[['span']] <- sister[['span']] - new_sister[['span']]
  node[['pd']] <- new_sister[['span']] + tip[['span']]
  node[['prdst']] <- sister[['prdst']] - new_sister[['span']]
  node[['prid']] <- sister[['prid']]
  node[['ptid']] <- node[['children']] <- c(tip[['id']], sister[['id']])
  tip[['pd']] <- 0
  tip[['prdst']] <- node[['prdst']] + tip[['span']]
  tip[['prid']] <- node[['id']]
  tree@nodelist[[tip[['id']]]] <- tip
  tree@nodelist[[node[['id']]]] <- node
  tree@nodelist[[new_sister[['id']]]] <- new_sister
  tree@nodelist[[new_parent[['id']]]] <- new_parent
  pres <- getNodePrid(tree, node[['id']])
  tree@nodelist[pres] <- lapply(tree@nodelist[pres],
                                    updatePre)
  .update(tree)
}
#TODO: add doc on pinning tips
pinTip <- function(tree, tip_id, lineage, end) {
  taxonyms <- unlist(lapply(tree@nodelist, function(n) n[['taxonym']]))
  for(i in length(lineage):1) {
    edges <- names(taxonyms)[which(taxonyms == lineage[i])]
    if(length(edges) == 0) {
      next
    }
    edges <- c(edges, unlist(sapply(edges, function(n) tree@nodelist[[n]][['ptid']])))
    edges <- edges[edges != tree@root]
    rngs <- getEdgesAge(tree, ids=edges)
    bool <- rngs[ ,'max'] > end
    if(any(bool)) {
      rngs <- rngs[bool, ]
      rngs[rngs[ ,'min'] <= end, "min"] <- end
      prbs <- rngs$max - rngs$min
      e <- as.vector(sample(rngs$edge, prob=prbs, size=1))
      e_i <- which(rngs$edge == e)
      start <- runif(min=rngs$min[e_i], max=rngs$max[e_i], n=1)
      if(i != length(lineage)) {
        tip_taxonym <- lineage[i+1]
      } else {
        tip_taxonym <- lineage[i]
      }
      tree <- addTip(tree, id=tip_id, sister=e, start=start, end=end,
                     tip_taxonym=tip_taxonym, parent_taxonym=lineage[i])
      break
    }
  }
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
