#' @name getNdAge
#' @title Get age
#' @description Return the age for \code{id}. Requires the known age of the tree to be provided.
#' @details Returns a numeric.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @param tree_age numeric value of known age of tree, tree['age'] if tree is up-to-date
#' @seealso
#' \code{\link{getNdsAge}}, 
#' \code{\link{getSpnAge}}, 
#' \code{\link{getSpnsAge}}, 
#' \code{\link{getPrnt}}, \code{\link{getTreeAge}}
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' # when did apes emerge?
#' # get parent id for all apes
#' prnt_id <- getPrnt(mammals, ids=c('Homo_sapiens', 'Hylobates_concolor'))
#' getNdAge(mammals, id=prnt_id, tree_age=mammals['age'])
getNdAge <- function(tree, id, tree_age) {
  prids <- .getSltPrids(tree@ndlst, FALSE)
  tree_age - .getNdPrdstsFrmLst(tree@ndlst, prids, id)
}

#' @name getSpnAge
#' @title Get age range
#' @description Return start and end ages for \code{id} from when it first appears to when it splits
#' @details Returns a dataframe.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @param tree_age numeric value of known age of tree, tree['age'] if tree is updated
#' @seealso
#' \code{\link{getNdAge}}, 
#' \code{\link{getNdsAge}}, 
#' \code{\link{getSpnsAge}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' data(mammals)
#' getSpnAge(mammals, id='Homo_sapiens', tree_age=mammals['age'])
getSpnAge <- function(tree, id, tree_age) {
  prids <- .getSltPrids(tree@ndlst, FALSE)
  end <- .getNdPrdstsFrmLst(tree@ndlst, prids, id)
  start <- end - tree@ndlst[[id]][['spn']]
  end <- tree_age - end
  start <- tree_age - start
  data.frame(spn=id, start, end)
}

#' @name getNdPrids
#' @title Get pre-nodes to root
#' @description Return node ids for connecting \code{id} to root.
#' @details Returns a vector. IDs are returned order from node ID to root.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsPrids}}, 
#' \code{\link{getNdPtids}}, 
#' \code{\link{getNdsPtids}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # get all nodes to root
#' getNdPrids(tree, id='t1')
getNdPrids <- function(tree, id) {
  prids <- .getSltPrids(tree@ndlst, FALSE)
  .getNdPridsFrmLst(tree@ndlst, prids, id)
}

#' @name getNdPtids
#' @title Get post-nodes to tips
#' @description Return node ids for connecting \code{id} to kids.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsPtids}}, 
#' \code{\link{getNdPrids}}, 
#' \code{\link{getNdsPrids}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # get all nodes from root to tip
#' getNdPtids(tree, id='n1')
# reduce dependence on the recursive, by getting prenodes
# tip ids to id
getNdPtids <- function(tree, id) {
  prids <- .getSltPrids(tree@ndlst, FALSE)
  .getNdPtidsFrmLst(tree@ndlst, prids, id)
}

#' @name getNdKids
#' @title Get children IDs
#' @description Return the node ids of all tips that descend from node.
#' @details Returns a vector
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsKids}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' # everyone descends from root
#' getNdKids(tree, id=tree['root'])
getNdKids <- function(tree, id) {
  prids <- .getSltPrids(tree@ndlst, FALSE)
  tids <- .getSltTids(tree@ndlst, FALSE)
  .getNdKidsFrmLst(tree@ndlst, prids, tids, id)
}

#' @name getNdPrdst
#' @title Get pre-distance
#' @description Return root to tip distance (prdst) for \code{id}
#' @details Sums the lengths of all branches from \code{id} to root.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsPrdst}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdPrdst(tree, id='t1')  # return the distance to root from t1
getNdPrdst <- function(tree, id) {
  prids <- .getSltPrids(tree@ndlst, FALSE)
  .getNdPrdstsFrmLst(tree@ndlst, prids, id)
}

#' @name getNdSlt
#' @title Get a node slot
#' @description Returns the value of named slot.
#' @details Returned object depends on name, either character, vector or numeric.
#' Default node slots are: id, spn, prid, ptid and txnym.
#' @param tree \code{TreeMan} object
#' @param slt_nm slot name
#' @param id node id
#' @seealso
#' \code{\link{getNdsSlt}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdSlt(tree, slt_nm='spn', id='t1')  # return span of t1
getNdSlt <- function(tree, slt_nm, id) {
  tree@ndlst[[id]][[slt_nm]]
}

#' @name getNdPD
#' @title Get phylogenetic diversity of node
#' @description Return summed value of all descending spns
#' @details Sums the lengths of all descending branches from a node.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsPD}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdPD(tree, id='n1')  # return PD of n1 which in this case is for the whole tree
getNdPD <- function(tree, id) {
  prids <- .getSltPrids(tree@ndlst, FALSE)
  .getNdPDFrmLst(tree@ndlst, prids, id)
}

#' @name getNdSstr
#' @title Get sister id
#' @description Returns the id of the sister(s) of node id given.
#' @details An error is raised if there is no sister (e.g. for the root).
#'  There can be more than one sister if tree is polytomous.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsSstr}}, 
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(10)
#' getNdSstr(tree, id='t1')
getNdSstr <- function(tree, id) {
  if(id == tree@root) {
    return(NULL)
  }
  .getNdSstrFrmLst(tree@ndlst, id)
}