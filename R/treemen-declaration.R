#' @name TreeMen-class
#' @title TreeMen-class
#' @aliases TreeMen-method
#' @param x \code{TreeMen} object
#' @param i tree index (integer or character)
#' @param object \code{TreeMen} object
#' @param max.level \code{str()} maximum level
#' @param ... additional tree objects
#' @param j missing
#' @param drop missing
#' @description S4 class for multiple phylogenetic trees
#' @slot treelst list of \code{TreeMan} objects
#' @slot ntips sum of tips per tree
#' @slot ntrees total number of trees
#' @exportClass TreeMen
#' @seealso 
#' \code{\link{cTrees}}
setClass('TreeMen', representation=representation(
  treelst='list',       # list of TreeMan objects
  ntips='numeric',       # sum of tips per tree
  ntrees='numeric'),     # number of trees in object
  validity=checkTreeMen)

# concatenate methods
.cMenToMen <- function(treemen_1, treemen_2) {
  treelst <- c(treemen_1@treelst, treemen_2@treelst)
  treemen_1@treelst <- treelst
  treemen_1@ntips <- treemen_1@ntips + treemen_2@ntips
  treemen_1@ntrees <- treemen_1@ntrees + treemen_2@ntrees
  treemen_1
}

.cMenToMan <- function(treemen, treeman) {
  treelst <- c(treemen@treelst, treeman)
  treemen@treelst <- treelst
  treemen@ntips <- treeman@ntips + treemen@ntips
  treemen@ntrees <- treemen@ntrees + 1
  treemen
}

.cMenToAny <- function(treemen, treeobj) {
  if(class(treeobj)[1] == "TreeMan") {
    treemen <- .cMenToMan(treemen, treeobj)
  } else if (class(treeobj)[1] == "TreeMen") {
    treemen <- .cMenToMen(treemen, treeobj)
  }
  treemen
}

.cTreeObjs <- function(treemen, treeobj, ...) {
  if(nargs() > 2) {
    treemen <- .cMenToAny(treemen, treeobj)
    treemen <- .cTreeObjs(treemen, ...)
  } else {
    treemen <- .cMenToAny(treemen, treeobj)
  }
  treemen
}

#' @title cTrees
#' @description Return \code{TreeMen} of concatenated trees.
#' @details Concatenate trees into single \code{TreeMen} object.
#' @param x \code{TreeMan} or \code{TreeMen} objects
#' @param ... more \code{TreeMan} or \code{TreeMen} objects
#' @seealso 
#' \code{\link{TreeMen-class}}, \code{\link{TreeMan-class}}, \code{\link{list-to-TreeMen}}
#' @examples 
#' library(treeman)
#' trees <- cTrees(randTree(10), randTree(10))
#' @export
setGeneric("cTrees", signature=c("x"),
           function(x, ...) {
             standardGeneric("cTrees")
           })
#' @rdname TreeMan-class
#' @exportMethod cTrees
setMethod("cTrees", c("TreeMan"),
          function(x, ...) {
            x <- list(x)
            x <- as(x, 'TreeMen')
            x <- .cTreeObjs(x, ...)
            x
            })
#' @rdname TreeMen-class
#' @exportMethod cTrees
setMethod('cTrees', c('TreeMen'),
          function(x, ...) {
            x <- .cTreeObjs(x, ...)
            x
          })

# Accessor methods
#' @rdname TreeMen-class
#' @exportMethod [[
setMethod('[[', c('TreeMen', 'ANY'),
          function(x, i) {
            if(!i %in% 1:length(x@treelst)) {
              stop(paste0(i, ' not in tree'))
            }
            x@treelst[[i]]
          })
#' @rdname TreeMen-class
#' @aliases TreeMen-method
#' Extract slots from a list of trees
#' @exportMethod [
setMethod('[', c('TreeMen', 'character', 'missing', 'missing'),
          function(x, i, j, ..., drop=TRUE) {
            if(!i %in% slotNames(x)) {
              stop(paste0(i, '  not in tree'))
            }
            slot(x, i)
          })

#' @name list-to-TreeMen
#' @title Convert list to a TreeMen
#' @description Return a \code{TreeMen} object from a list of \code{TreeMans}
#' @seealso 
#' \code{\link{TreeMen-class}}
#' @examples 
#' library(treeman)
#' trees <- list('tree_1'=randTree(10), 'tree_2'=randTree(10))
#' trees <- as(trees, 'TreeMen')
# Conversion method
setAs(from="list", to="TreeMen", def=function(from, to) {
  ntips <- sum(unlist(lapply(from, function(tree) tree@ntips)))
  ntrees <- length(from)
  new(to, treelst=from, ntips=ntips, ntrees=ntrees)
  })

# display methods
#' @rdname TreeMen-class
#' @exportMethod as.character
setMethod('as.character', c('x'='TreeMen'),
          function(x) {
            paste0('TreeMen Object of [', x@ntrees,'] trees')
          })
#' @rdname TreeMen-class
#' @exportMethod show
setMethod('show', 'TreeMen',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname TreeMen-class
#' @exportMethod str
setMethod('str', c('object'='TreeMen'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
#' @rdname TreeMen-class
#' @exportMethod print
setMethod('print', c('x'='TreeMen'),
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname TreeMen-class
#' @exportMethod summary
setMethod('summary', c('object'='TreeMen'),
          function(object){
            msg <- 'Trees (TreeMen Object):\n'
            msg <- paste0(msg, '  + ', object@ntrees, ' trees\n')
            msg <- paste0(msg, '  + ', object@ntips, ' tips\n')
            cat(msg)
          })