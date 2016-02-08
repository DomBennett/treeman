# TODO: write Newick

# TODO: readTree doc
readTree <- function(file=NULL, text=NULL, ...) {
  if(!is.null(file)) {
    trstr <- scan(file, what="raw", quiet=TRUE)
  } else {
    trstr <- text
  }
  if(length(trstr) > 1) {
    trees <- mlply(trstr, .fun=.readTree, ...)
    tree <- as(trees, 'TreeMen')
  } else {
    tree <- .readTree(trstr)
  }
  tree
}

.readTree <- function(trstr) {
  cuts <- gregexpr("(\\(|\\)|,|;)", trstr)[[1]]
  cuts <- c(cuts[1], cuts[2:length(cuts)] - cuts[1:(length(cuts)-1)])
  rdrenv <- .getRdrEnv(trstr)
  l_data <- data.frame(end_pos=cuts, stringsAsFactors=FALSE)
  m_ply(l_data, .mkNdLst, rdrenv=rdrenv)
  .addRoot(rdrenv)
  if(length(rdrenv[['root']]) == 0) {
    ndlst <- .globalUpdateKids(rdrenv[['nodelist']])
  } else {
    ndlst <- .globalUpdateAll(rdrenv[['nodelist']])
  }
  tree <- new('TreeMan', nodelist=ndlst, root=rdrenv[['root']])
  .updateSlots(tree)
}

# set-up reader env
.getRdrEnv <- function(trstr) {
  rdrenv <- new.env()
  rdrenv$wspn <- grepl(':', trstr)
  rdrenv$trstr <- trstr
  rdrenv$nodelist <- list()
  rdrenv$prnds <- list()
  rdrenv$cntr <- 0L
  rdrenv$i <- 0L
  rdrenv$nxt_is_intrnl <- FALSE
  rdrenv
}

# extract ID and span from ndstr
.getIDandSpan <- function(ndstr, nints, wspn) {
  nd <- .mkNd(id='', wspn)
  ndstr <- gsub("(\\(|\\)|\\;|,)", "", ndstr)
  ndstr <- strsplit(ndstr, ":")[[1]]
  if(length(ndstr) > 1) {
    nd[['span']] <- as.numeric(ndstr[2])
  }
  if(length(ndstr) == 0 || ndstr[1] == "") {
    nd[['id']] <- paste0("n", nints)
  } else {
    nd[['id']] <- ndstr[1]
  }
  nd
}

.mkNd <- function(id, wspn) {
  nd <- list()
  nd[['id']] <- id
  if(wspn) {
    nd[['span']] <- nd[['prdst']] <- nd[['pd']] <- 0
  }
  nd
}

.getPrid <- function(prnds) {
  unlist(lapply(prnds[length(prnds):1], function(nd) nd[['id']]))
}

# cut trstr and generate nodelist
.mkNdLst <- function(end_pos, rdrenv) {
  ndstr <- substr(rdrenv$trstr, 1, end_pos)
  if(grepl("^\\(", ndstr)) {
    rdrenv$cntr <- rdrenv$cntr + 1
    rdrenv$i <- rdrenv$i + 1
    rdrenv$prnds[[rdrenv$i]] <- .mkNd(id=paste0("n", rdrenv$cntr),
                                      wspn=rdrenv$wspn)
  } else {
    nd <- .getIDandSpan(ndstr, rdrenv$cntr, rdrenv$wspn)
    if(rdrenv$nxt_is_intrnl) {
      # TODO: utilise nd$id, e.g. node labels are taxonym or support
      nd$id <- rdrenv$prnds[[rdrenv$i]][['id']]
      nd$ptid <- rdrenv$prnds[[rdrenv$i]][['ptid']]
      rdrenv$prnds <- rdrenv$prnds[-rdrenv$i]
      rdrenv$i <- rdrenv$i - 1
      nd$prid <- .getPrid(rdrenv$prnds)
      rdrenv$nodelist[[nd[['id']]]] <- nd
      rdrenv$nxt_is_intrnl <- FALSE
    } else {
      nd[['prid']] <- .getPrid(rdrenv$prnds)
      rdrenv$nodelist[[nd[['id']]]] <- nd
    }
    if(length(rdrenv$prnds) > 0) {
      rdrenv$prnds[[rdrenv$i]][['ptid']] <-
        c(rdrenv$prnds[[rdrenv$i]][['ptid']], nd[['id']])
    }
    if(grepl("\\)$", ndstr)) {
      rdrenv$nxt_is_intrnl <- TRUE
    }
  }
  rdrenv$trstr <- substr(rdrenv$trstr, end_pos+1, nchar(rdrenv$trstr))
  NULL
}

.addRoot <- function(rdrenv) {
  root_i <- which(unlist(lapply(rdrenv$nodelist, function(n) length(n[['prid']]) == 0)))
  if(length(root_i) > 0) {
    rdrenv$nodelist[[root_i]][['prid']] <- NULL
    rdrenv$nodelist[[root_i]][['span']] <- 0
    rdrenv$root <- names(rdrenv$nodelist)[root_i]
  } else {
    rdrenv$root <- character()
  }
  NULL
}