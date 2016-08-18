#' phylo class
#'
#' @name phylo-class
#' @aliases phylo
#'
#' @exportClass phylo
setOldClass('phylo')

#' @name TreeMan-to-phylo
#' @title Convert TreeMan to phylo
#' @description Return ape's \code{phylo} from a \code{TreeMan}
#' @seealso 
#' \code{\link{phylo-to-TreeMan}},
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

# TODO: multiphylo conversions