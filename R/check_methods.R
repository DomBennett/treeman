#' @name fastCheckTreeMan
#' @title Check if tree is correct, fast!
#' @description Return T/F if tree is a true \code{TreeMan} object
#' @details Whenever a tree is first initiated this check is used.
#' For more detailed checking use \code{checkNdlst}.
#' @param object \code{TreeMan} object
#' @seealso
#' \code{\link{checkNdlst}}, \code{\link{checkTreeMen}}
#' @export
fastCheckTreeMan <- function(object) {
  kwn_ids <- names(object@ndlst)
  ids <- unlist(sapply(object@ndlst, function(x) x[c('id', 'ptid', 'prid')]))
  test <- all(ids %in% kwn_ids) & object@root %in% kwn_ids
  # check hierarchy through prinds
  prinds <- object@prinds
  if(length(prinds) > 0) {
    # only root node should equal its index
    prind_test <- sum(prinds == 1:length(prinds)) == 1
    # all internal nodes should occur more than once (twice for bifurcating trees)
    prind_test <- all(table(prinds) > 1) & prind_test
    test <- test & prind_test
  }
  test
}

#' @name checkNdlst
#' @title Check if ndlst is correct
#' @description Return T/F fpr \code{ndlst} consistency
#' @details Tests whether each node in tree points to valid other node IDs. Also
#' ensures `spn` and `root` are correct. Reports nodes that have errors.
#' @param ndlst \code{ndlst}
#' @param root root ID
#' @seealso
#' \code{\link{fastCheckTreeMan}}, \code{\link{checkTreeMen}}
#' @export
#' @examples
#' library(treeman)
#' tree <- randTree(100)
#' (checkNdlst(tree@ndlst, tree@root))
checkNdlst <- function(ndlst, root) {
  .check <- function(nd) {
    # must have id
    test_id <- is.character(nd[['id']]) & 'id' %in% names(nd)
    # id must contain no special characters
    test_spcl_chrs <- test_id && !grepl('[^0-9a-zA-Z_]', nd[['id']])
    # txnyms
    test_txnym <- TRUE
    if('txnym' %in% names(nd)) {
      test_txnym <- is.character(nd[['txnym']])
      for(txnym in nd[['txnym']]) {
        test_txnym <- test_txnym && !grepl('[^0-9a-zA-Z_]', txnym)
      }
    }
    # must have either prid/ptid or both
    test_slts <- ('ptid' %in% names(nd) | 'prid' %in% names(nd))
    test_valid_nd <- nd[['id']] %in% nds  # nd id must be known
    # prid and ptids must be known
    test_prid <- is.character(nd[['prid']]) & nd[['prid']] %in% nds
    test_ptid <- is.character(nd[['ptid']]) > 0 & all(nd[['ptid']] %in% nds)
    # spns must be 0 or more
    test_spn <- is.numeric(nd[['spn']]) && nd[['spn']] >= 0
    # test self-reference
    test_sr <- !nd[['prid']] %in% nd['ptid']
    # test root is never a ptid, proxy for circularity
    test_circ <- !rid %in% nd[['ptid']]
    # only root is self-referential
    test_root <- rid != nd[['id']] |
      (rid == nd[['id']] & rid == nd[['prid']])
    bool <- test_id & test_valid_nd & test_prid &
      test_ptid & test_sr & test_circ & test_slts &
      test_spcl_chrs & test_txnym
    if(length(bool) > 0 && bool) {
      return(TRUE)
    }
    FALSE
  }
  nds <- names(ndlst)
  rid <- root
  nd_checks <- sapply(ndlst, .check)
  if(!all(nd_checks)) {
    msg <- 'These nodes are invalid:\n'
    bad <- which(!nd_checks)
    for(i in bad[-length(bad)]) {
      msg <- paste0(msg, nds[i], ', ')
    }
    msg <- paste0(msg, nds[bad[length(bad)]], '\n')
    cat(msg, '\n')
    return(FALSE)
  }
  if(!rid %in% nds) {
    msg <- paste0("Root node `", rid, '` not in ndlst\n')
    cat(msg, '\n')
    return(FALSE)
  }
  TRUE
}

#' @name checkTreeMen
#' @title Check if trees are correct
#' @description Return T/F if trees is a true \code{TreeMen} object
#' @details Tests whether all trees in object are \code{TreeMan} objects
#' @param object \code{TreeMen} object
#' @seealso
#' \code{\link{checkNdlst}}
#' @export
checkTreeMen <- function(object) {
  .check <- function(i, invlds) {
    if(class(object@treelst[[i]])[1] != "TreeMan") {
      invlds <<- c(i, invlds)
    }
    NULL
  }
  invlds <- NULL
  sapply(1:length(object@treelst), .check, invlds=invlds)
  if(length(invlds) > 0) {
    for(i in invlds) {
      cat("[", i, "] in treelst is not a TreeMan object\n", sep="")
    }
    return(FALSE)
  }
  TRUE
}