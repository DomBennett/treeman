#' @name fastCheckTreeMan
#' @title Check if tree is correct, fast!
#' @description Return T/F if tree is a true \code{TreeMan} object
#' @details Whenever a tree is first initiated this check is used.
#' For more detailed checking use \code{checkTreeMan}.
#' @param object \code{TreeMan} object
#' @seealso
#' \code{\link{checkTreeMan}}, \code{\link{checkTreeMen}}
#' @export
fastCheckTreeMan <- function(object) {
  kwn_ids <- names(object@ndlst)
  ids <- unlist(sapply(object@ndlst, function(x) x[c('id', 'ptid', 'prid')]))
  test <- all(ids %in% kwn_ids) & object@root %in% kwn_ids
  test
}

#' @name checkTreeMan
#' @title Check if tree is correct
#' @description Return T/F if tree is a true \code{TreeMan} object
#' @details Tests whether each node in tree points to valid other node IDs. Also
#' ensures `spn` and `root` are correct. Reports nodes that have errors.
#' @param object \code{TreeMan} object
#' @seealso
#' \code{\link{fastCheckTreeMan}}, \code{\link{checkTreeMen}}
#' @export
checkTreeMan <- function(object) {
  # TODO: use prids as vector to test for circularity
  .check <- function(nd) {
    test_id <- is.character(nd[['id']]) & 'id' %in% names(nd)  # must have id
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
      test_ptid & test_sr & test_circ & test_slts
    if(length(bool) > 0 && bool) {
      return(TRUE)
    }
    FALSE
  }
  nds <- names(object@ndlst)
  rid <- object@root
  nd_checks <- sapply(object@ndlst, .check)
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
    msg <- paste0("Root node `", rid, '` not in tree.\n')
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
#' \code{\link{checkTreeMan}}
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