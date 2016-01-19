# TODO: write Newick

# TODO: readTree doc
readTree <- function(file=NULL, tree_string=NULL) {
  if(!is.null(file)) {
    trstr <- scan(file, what="raw", quiet=TRUE)
  } else {
    trstr <- tree_string
  }
  #trstr <- sub(";", "", trstr)
  cuts <- gregexpr("(\\(|\\)|,|;)", trstr)[[1]]
  cuts <- c(cuts[1], cuts[2:length(cuts)] - cuts[1:(length(cuts)-1)])
  rdrenv <- .getRdrEnv(trstr)
  sapply(cuts, .mkNdLst, rdrenv=rdrenv)
  # TODO merge these into a single function to make faster
  .addRoot(rdrenv)
  .addChildren(rdrenv)
  .addPredist(rdrenv)
  .addPD(rdrenv)
  tree <- new('TreeMan', nodelist=rdrenv$nodelist, root=rdrenv$root)
  .update(tree)
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
  root_i <- which(unlist(lapply(rdrenv$nodelist, function(n) n$prenode == "root")))
  if(length(root_i) > 0) {
    rdrenv$nodelist[[root_i]]$prenode <- NULL
    rdrenv$nodelist[[root_i]]$span <- 0
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
    if(!is.null(nd$prenode)) {
      d <- nd$span + d
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
  rdrenv$prenodes[[1]] <- list()
  rdrenv$prenodes[[1]]$id <- "root"
  rdrenv$cntr <- 0L
  rdrenv$i <- 1L
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
    nd$span <- NULL
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
    rdrenv$cntr <- rdrenv$cntr + 1
    rdrenv$i <- rdrenv$i + 1
    rdrenv$prenodes[[rdrenv$i]] <- list()
    rdrenv$prenodes[[rdrenv$i]]$id <- paste0("n", rdrenv$cntr)
  } else {
    nd <- .getIDandSpan(ndstr, rdrenv$cntr)
    if(rdrenv$nxt_is_intrnl) {
      # TODO: utilise nd$id, e.g. node labels are taxonym or support
      nd$id <- rdrenv$prenodes[[rdrenv$i]]$id
      rdrenv$prenodes[[rdrenv$i]]$span <- nd$span
      rdrenv$prenodes[[rdrenv$i]]$prenode <-
        rdrenv$prenodes[[rdrenv$i-1]]$id
      rdrenv$nodelist[[nd$id]] <- rdrenv$prenodes[[rdrenv$i]]
      rdrenv$prenodes <- rdrenv$prenodes[-rdrenv$i]
      rdrenv$i <- rdrenv$i - 1
      rdrenv$nxt_is_intrnl <- FALSE
    } else {
      nd$prenode <- rdrenv$prenodes[[rdrenv$i]]$id
      rdrenv$nodelist[[nd$id]] <- nd
    }
    if(length(rdrenv$prenodes) > 0) {
      rdrenv$prenodes[[rdrenv$i]]$postnode <-
        c(rdrenv$prenodes[[rdrenv$i]]$postnode, nd$id)
    }
    if(grepl("\\)$", ndstr)) {
      rdrenv$nxt_is_intrnl <- TRUE
    }
  }
  rdrenv$trstr <- substr(rdrenv$trstr, end_pos+1, nchar(rdrenv$trstr))
  NULL
}
