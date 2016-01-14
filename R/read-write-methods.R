# TODO: write Newick

# TODO: readTree doc
readTree <- function(file=NULL, tree_string=NULL) {
  if(!is.null(file)) {
    trstr <- scan(file, what="raw")
  } else {
    trstr <- tree_string
  }
  cuts <- gregexpr("(\\(|\\)|,|;)", trstr)[[1]]
  cuts <- c(cuts[1], cuts[2:length(cuts)] - cuts[1:(length(cuts)-1)])
  rdrenv <- .getRdrEnv(trstr)
  sapply(cuts, .mkNdLst, rdrenv=rdrenv)
  # TODO handle trees without branch lengths
  .addRoot(rdrenv)
  .addChildren(rdrenv)
  .addPredist(rdrenv)
  .addPD(rdrenv)
  tree <- new ('TreeMan', nodelist=rdrenv$nodelist, root=rdrenv$root)
  tree <- .update(tree)
  tree
}

# addPD
.addPD <- function(rdrenv) {
  assgn <- function(i) {
    rdrenv$nodelist[[i]]$pd <- 0
    NULL
  }
  add <- function(prndid, rdrenv, ndid=prndid) {
    if(!is.null(rdrenv$nodelist[[prndid]]$postnode)) {
      rdrenv$nodelist[[prndid]]$pd <- rdrenv$nodelist[[prndid]]$pd +
        rdrenv$nodelist[[ndid]]$span
    }
    if(!is.null(rdrenv$nodelist[[prndid]]$prenode)) {
      add(rdrenv$nodelist[[prndid]]$prenode, rdrenv, ndid)
    }
    NULL
  }
  sapply(1:length(rdrenv$nodelist), assgn)
  sapply(names(rdrenv$nodelist), add, rdrenv=rdrenv)
  NULL
}

# add children
.addChildren <- function(rdrenv) {
  add <- function(ndid, rdrenv, tpid=ndid) {
    if(!is.null(rdrenv$nodelist[[ndid]]$postnode)) {
      rdrenv$nodelist[[ndid]]$children <- c(tpid, rdrenv$nodelist[[ndid]]$children)
    }
    if(!is.null(rdrenv$nodelist[[ndid]]$prenode)) {
      add(rdrenv$nodelist[[ndid]]$prenode, rdrenv, tpid)
    }
    NULL
  }
  tips <- sapply(rdrenv$nodelist, function(n) length(n$postnode) == 0)
  tips <- names(tips)[tips]
  sapply(tips, add, rdrenv=rdrenv)
  NULL
}

.addRoot <- function(rdrenv) {
  root_i <- which(unlist(lapply(rdrenv$nodelist, function(n) n$prenode == "n0")))
  if(length(root_i) > 0) {
    rdrenv$nodelist[[root_i]]$prenode <- NULL
    rdrenv$root <- names(rdrenv$nodelist)[root_i]
  } else {
    rdrenv$root <- character()
  }
  NULL
}

# add predists
.addPredist <- function(rdrenv) {
  calc <- function(nd, d) {
    nd <- rdrenv$nodelist[[nd]]
    d <- nd$span + d
    if(!is.null(nd$prenode)) {
      d <- calc(nd$prenode, d)
    }
    d
  }
  assgn <- function(i) {
    rdrenv$nodelist[[i]]$predist <- ds[[i]]
    NULL
  }
  ds <- sapply(names(rdrenv$nodelist), calc, d=0)
  sapply(1:length(ds), assgn)
}

# set-up reader env
.getRdrEnv <- function(trstr) {
  rdrenv <- new.env()
  rdrenv$trstr <- trstr
  rdrenv$nodelist <- list()
  rdrenv$prenodes <- list()
  rdrenv$i <- 0L
  rdrenv$nxt_is_intrnl <- FALSE
  rdrenv
}

# extract ID and span from ndstr
.getIDandSpan <- function(ndstr, nints) {
  nd <- list()
  ndstr <- gsub("(\\(|\\)|\\;|,)", "", ndstr)
  ndstr <- strsplit(ndstr, ":")[[1]]
  if(length(ndstr) > 1) {
    nd$span <- as.numeric(ndstr[2])
  } else {
    nd$span <- 0
  }
  if(length(ndstr) == 0 || ndstr[1] == "") {
    nd$id <- paste0("n", nints)
  } else {
    nd$id <- ndstr[1]
  }
  nd
}

# cut trstr and generate nodelist
.mkNdLst <- function(end_pos, rdrenv) {
  ndstr <- substr(rdrenv$trstr, 1, end_pos)
  if(grepl("^\\(", ndstr)) {
    rdrenv$i <- length(rdrenv$prenodes) + 1
    rdrenv$prenodes[[rdrenv$i]] <- list()
  } else {
    nd <- .getIDandSpan(ndstr, rdrenv$i)
    if(rdrenv$nxt_is_intrnl) {
      # TODO: utilise nd$id
      nd$id <- paste0('n', rdrenv$i)
      rdrenv$prenodes[[rdrenv$i]]$id <- nd$id
      rdrenv$prenodes[[rdrenv$i]]$span <- nd$span
      rdrenv$prenodes[[rdrenv$i]]$prenode <- paste0('n', rdrenv$i - 1)
      rdrenv$nodelist[[nd$id]] <- rdrenv$prenodes[[rdrenv$i]]
      rdrenv$i <- rdrenv$i - 1
      rdrenv$nxt_is_intrnl <- FALSE
    } else {
      nd$prenode <- paste0('n', rdrenv$i)
      rdrenv$nodelist[[nd$id]] <- nd
    }
    if(rdrenv$i > 0) {
      rdrenv$prenodes[[rdrenv$i]]$postnode <-
        c(rdrenv$prenodes[[rdrenv$i]]$postnode, nd$id)
    }
    if(grepl("(\\)|;)$", ndstr)) {
      rdrenv$nxt_is_intrnl <- TRUE
    }
  }
  rdrenv$trstr <- substr(rdrenv$trstr, end_pos+1, nchar(rdrenv$trstr))
}
