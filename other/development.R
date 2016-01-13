# TODO:
# -- read/write newick


library(treeman)
readTree <- function(file=NULL, tree_string=NULL) {
  getIDandSpan <- function(ndstr, n_intnds) {
    nd <- list()
    ndstr <- gsub("(\\(|\\)|\\;|,)", "", ndstr)
    ndstr <- strsplit(ndstr, ":")[[1]]
    # TODO: calculate predist
#     if(length(ndstr) > 1) {
#       nd$span <- as.numeric(ndstr[2])
#     } else {
#       nd$span <- 0
#     }
    if(length(ndstr) == 0 || ndstr[1] == "") {
      nd$id <- paste0("n", n_intnds)
    } else {
      nd$id <- ndstr[1]
    }
    nd
  }
  mkNdLst <- function(end_pos) {
    ndstr <- substr(trstr, 1, end_pos)
    if(grepl("^\\(", ndstr)) {
      i <- length(prenodes) + 1
      prenodes[[i]] <- list()
    } else {
      nd <- getIDandSpan(ndstr, i)
      if(nxt_is_intrnl) {
        # TODO: utilise nd$id
        nd$id <- paste0('n', i)
        prenodes[[i]]$id <- nd$id
        prenodes[[i]]$span <- nd$span
        prenodes[[i]]$prenode <- paste0('n', i - 1)
        nodelist[[nd$id]] <- prenodes[[i]]
        i <- i - 1
        nxt_is_intrnl <- FALSE
      } else {
        nd$prenode <- paste0('n', i)
        nodelist[[nd$id]] <- nd
      }
      if(i > 0) {
        prenodes[[i]]$postnode <- c(prenodes[[i]]$postnode,
                                    nd$id)
      }
      if(grepl("(\\)|;)$", ndstr)) {
        nxt_is_intrnl <- TRUE
      }
    }
    trstr <<- substr(trstr, end_pos+1, nchar(trstr))
    i <<- i
    print(i)
    nxt_is_intrnl <<- FALSE
    print(nxt_is_intrnl)
    prenodes <<- prenodes
    nodelist <<- nodelist
  }
  if(!is.null(file)) {
    trstr <- scan(file, what="raw")
  } else {
    trstr <- tree_string
  }
  nodelist <- list()
  prenodes <- list()
  i <- 0
  nxt_is_intrnl <- FALSE
  cuts <- gregexpr("(\\(|\\)|,|;)", trstr)[[1]]
  cuts <- c(cuts[1], cuts[2:length(cuts)] - cuts[1:(length(cuts)-1)])
  sapply(cuts, mkNdLst)
  root_i <- which(unlist(lapply(nodelist, function(n) n$prenode == "n0")))
  if(length(root_i) > 0) {
    nodelist[[root_i]]$prenode <- NULL
    root <- nodelist[[root_i]]$id
    tree <- new ('TreeMan', nodelist=nodelist, root=root)
  } else {
    tree <- new ('TreeMan', nodelist=nodelist)
  }
  tree <- .update(tree)
}

readTree(tree_string="(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5);")
