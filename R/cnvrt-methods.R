#' phylo class
#'
#' @name phylo-class
#' @aliases phylo
#'
#' @exportClass phylo
setOldClass('phylo')

#' multiPhylo class
#'
#' @name multiPhylo-class
#' @aliases multiPhylo
#'
#' @exportClass multiPhylo
setOldClass('multiPhylo')

#' @name TreeMan-to-phylo
#' @title Convert TreeMan to phylo
#' @description Return ape's \code{phylo} from a \code{TreeMan}
#' @seealso 
#' \code{\link{phylo-to-TreeMan}},
#' \code{\link{TreeMen-to-multiPhylo}}
#' \code{\link{multiPhylo-to-TreeMen}}
#' \code{\link{TreeMan-class}}
#' @examples 
#' library(treeman)
#' library(ape)
#' tree <- randTree(10)
#' tree <- as(tree, 'phylo')
setAs(from="phylo", to="TreeMan", def=function(from, to) {
  ape::write.tree(from, file='temp.tre')
  tree <- readTree(file='temp.tre')
  file.remove('temp.tre')
  return(tree)
})

#' @name phylo-to-TreeMan
#' @title Convert phylo to TreeMan
#' @description Return a \code{TreeMan} from ape's \code{phylo}
#' @seealso
#' \code{\link{TreeMan-to-phylo}},
#' \code{\link{TreeMen-to-multiPhylo}}
#' \code{\link{multiPhylo-to-TreeMen}}
#' \code{\link{TreeMan-class}}
#' @examples 
#' library(treeman)
#' library(ape)
#' tree <- compute.brlen(rtree(10))
#' tree <- as(tree, 'TreeMan')
setAs(from="TreeMan", to="phylo", def=function(from, to) {
  writeTree(from, file='temp.tre')
  tree <- ape::read.tree(file='temp.tre')
  file.remove('temp.tre')
  return(tree)
})

#' @name multiPhylo-to-TreeMen
#' @title Convert multiPhylo to TreeMen
#' @description Return a \code{TreeMen} from ape's \code{mutlPhylo}
#' @seealso
#' \code{\link{TreeMan-to-phylo}},
#' \code{\link{phylo-to-TreeMan}},
#' \code{\link{TreeMen-to-multiPhylo}}
#' \code{\link{TreeMan-class}}
#' @examples 
#' library(treeman)
#' library(ape)
#' trees <- c(rtree(10), rtree(10), rtree(10))
#' trees <- as(trees, 'TreeMen')
setAs(from="multiPhylo", to="TreeMen", def=function(from, to) {
  ape::write.tree(from, file='temp.tre')
  tree <- readTree(file='temp.tre')
  file.remove('temp.tre')
  return(tree)
})

#' @name TreeMen-to-multiPhylo
#' @title Convert TreeMen to multiPhylo
#' @description Return ape's \code{multiPhylo} from a \code{TreeMen}
#' @seealso 
#' \code{\link{TreeMan-to-phylo}},
#' \code{\link{phylo-to-TreeMan}},
#' \code{\link{multiPhylo-to-TreeMen}}
#' \code{\link{TreeMan-class}}
#' @examples 
#' library(treeman)
#' library(ape)
#' trees <- cTrees(randTree(10), randTree(10), randTree(10))
#' trees <- as(trees, 'multiPhylo')
setAs(from="TreeMen", to="multiPhylo", def=function(from, to) {
  writeTree(from, file='temp.tre')
  tree <- ape::read.tree(file='temp.tre')
  file.remove('temp.tre')
  return(tree)
})