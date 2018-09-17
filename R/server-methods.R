#' @name searchTxnyms
#' @title Get node labels based on online taxonomic database
#' @description Return names of each node in tree based on searching tip labels
#'  through Global Names Resolver \url{http://resolver.globalnames.org/} in NCBI.
#' @details For each node, all the descendants are searched, the taxonomic lineages returned and
#' then searched to find the lowest shared name.
#' All the tip labels are searched against a specified taxonomic database through the GNR and NCBI.
#' (So far only tested with NCBI database.)
#' Use the infer argument to ensure a taxonym is returned for all nodes. If infer is true,
#' all nodes without an identifed taxonym adopt the taxonym of their parent.
#' Will raise a warning if connection fails and will return NULL.
#' @param tree TreeMan object
#' @param cache T/F, create a local cache of downloaded names?
#' @param parent specify parent of all names to prevent false names
#' @param clean T/F, ensure returned names contain no special characters?
#' @param infer T/F, infer taxonyms for unfound nodes?
#' @seealso
#' \code{\link{taxaResolve}}, \code{\link{setTxnyms}}, \code{\link{getNdsFrmTxnyms}}
#' @export
#' @examples
#' tree <- randTree(8)
#' new_tids <- c("Gallus_gallus", "Aileuropoda_melanoleucha", "Ailurus_fulgens",
#' "Rattus_rattus", "Mus_musculus", "Gorilla_gorilla", "Pan_trogoldytes", "Homo_sapiens")
#' tree <- setNdsID(tree, tree['tips'], new_tids)
#' nd_labels <- searchTxnyms(tree)
#' print(nd_labels)
# TODO: add compatibility with other GNR datasources
# TODO: catalogue of life, unlike NCBI, does not keep lineages and rank lengths constant between names
searchTxnyms <- function (tree, cache=FALSE, parent=NULL, clean=TRUE,
                          infer=TRUE) {
  # Use GNR to label all nodes in a phylogeny
  # first replace all _ with spaces
  tip_labels <- gsub ('_', ' ', tree@tips)
  nids <- tree@nds
  nd_labels <- rep(NA, tree@nall)
  names(nd_labels) <- tree@all
  taxa_res <- taxaResolve(tip_labels, datasource=4, cache=cache, parent=parent)
  if(is.null(taxa_res)) {
    return(NULL)
  }
  # for tips use the first word of the name
  nd_labels[tree@tips] <- vapply(strsplit(tip_labels, "\\s+"), function(x) x[1],
                                 character(1))
  nds_kids <- getNdsKids(tree, ids=nids)
  for(nid in nids) {
    kids <- gsub("_" , " ", nds_kids[[nid]])
    genus_names <- vapply(strsplit(kids, "\\s+"), function(x) x[1],
                          character(1))
    if(all(genus_names == genus_names[1])) {
      nd_labels[nid] <- genus_names[1]
    } else {
      lineages <- as.character(taxa_res[taxa_res$search.name %in% kids, "lineage"])
      lineages <- strsplit (lineages, "\\|")
      lineages <- lineages[!is.na (lineages)]
      if(length (lineages) > 1) {
        nd_labels[nid] <- .findClade(lineages)
      }
    }
  }
  if(clean) {
    nd_labels <- gsub('\\s', '_', nd_labels)
    nd_labels <- gsub('[^a-zA-Z_0-9]', '', nd_labels)
  }
  if(infer & any(is.na(nd_labels))) {
    if(sum(is.na(nd_labels))/length(nd_labels) > 0.5) {
      message('Fewer than 50% of nodes identified,',
              ' not attempting inference of remaning nodes.')
    } else {
      for(i in which(is.na(nd_labels))) {
        prids <- getNdPrids(tree, names(nd_labels)[i])
        pssbls <- nd_labels[prids]
        pssbls <- pssbls[!is.na(pssbls)]
        if(length(pssbls) > 0) {
          nd_labels[[i]] <- pssbls[[1]]
        }
      }
    }
  }
  nd_labels
}

#' @name taxaResolve
#' @title Resolve taxonmic names online
#' @description Resolve taxonomic names via the Global Names Resolver.
#' @details Returns dataframe containing GNR metadata for each name wames
#' that cannot be resolved are returned as NA. Various datasources are 
#' available, see \url{http://resolver.globalnames.org/data_sources} for a
#' list and IDs. Default is 4 for NCBI.
#' Will raise a warning if connection fails and will return NULL.
#' @param nms vector of names
#' @param batch size of the batches to be queried
#' @param datasource ID number of the datasource
#' @param genus boolean, if true will search against GNR with just the genus
#'  name for names that failed to resolve using the full species name
#' @param cache T/F, create a local cache of downloaded names?
#' @param parent specify parent of all names to prevent false names
#' @seealso
#' \code{\link{searchTxnyms}}, \code{\link{setTxnyms}}, \code{\link{getNdsFrmTxnyms}}
#' @export
#' @examples
#' my_lovely_names <- c ('Gallus gallus', 'Pongo pingu', 'Homo sapiens',
#'                       'Arabidopsis thaliana', 'Macaca thibetana', 'Bacillus subtilis')
#' res <- taxaResolve (nms=my_lovely_names)
#' length(colnames(res))  # 10 different metadata for returned names including original search name
#' # let's look at the lineages
#' lineages <- strsplit(as.vector(res$lineage), '\\|')
#' print(lineages[[6]])  # the bacteria has far fewer taxonomic levels
# NOTE. Originally built for MTT
taxaResolve <- function (nms, batch=100, datasource=4, genus=TRUE,
                         cache=FALSE, parent=NULL) {
  .replace <- function (i, slot.name) {
    # controlled extraction
    element <- data[[i]]$result[[1]][[slot.name]]
    if (!is.null (element)) {
      res <- element
    } else {
      res <- NA
    }
    res
  }
  batchResolve <- function (batch.nms) {
    #create query from nms
    url <- "http://resolver.globalnames.org/name_resolvers.json?"
    data_source_ids <- paste0 ("&data_source_ids=", datasource)
    nms2 <- paste0 ("names=", paste0 (stringr::str_replace_all (
      batch.nms, " ", "+"), collapse = "|"))
    query <- paste (plyr::compact (list (url, nms2,
                                         data_source_ids)), collapse = "")
    #search via API
    data <- .safeFromJSON(query)$data
    return(data)
  }
  # avoid names -- names exist on database but mean nothing
  avoid <- c ('unidentified')
  # make sure names don't have '_'
  trms <- gsub ('_', ' ', nms)
  # remove all non alphanumerics
  trms <- gsub ('\\W+', ' ', trms)
  # remove trailing whitespace
  trms <- gsub ('^\\s+|\\s+$', '', trms)
  # any missing trms, replace with stubs
  trms[trms==''] <- 'invalid'
  deja_vues <- rep(FALSE, length(trms))
  data <- vector("list", length=length(trms))
  names(data) <- nms
  if(cache) {
    if(!file.exists("gnr_cache")) {
      dir.create("gnr_cache")
    }
    for(i in 1:length(nms)) {
      fp <- file.path("gnr_cache", paste0(nms[i], ".RData"))
      if(file.exists(fp)) {
        load(fp)
        data[[nms[i]]] <- nd
        deja_vues[i] <- TRUE
      }
    }
  }
  if(sum(!deja_vues) > 0) {
    # Split nms into batch sized chunks
    #  http://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
    x <- seq_along (trms[!deja_vues])
    btrms <- split (trms[!deja_vues], ceiling (x/batch))
    bnms <- split (nms[!deja_vues], ceiling (x/batch))
    for (i in 1:length(btrms)) {
      temp.data <- batchResolve(btrms[[i]])
      if (is.null(temp.data)) {
        return(NULL)
      }
      data[bnms[[i]]] <- temp.data
    }
  }
  #transform results into output
  search.name <- name.string <- canonical.form <-
    lineage <- lineage.ids <- rank <- taxid <-
    match.type <- prescore <- score <- rep (NA, length (nms))
  for (i in 1:length (data)){
    parent_test <- TRUE
    nd <- data[[i]]
    if(cache & !deja_vues[i]) {
      fp <- file.path("gnr_cache", paste0(names(data)[i], ".RData"))
      save(nd, file=fp)
    }
    if (!'results' %in% names (nd)){
      search.name[i] <- nms[i]
    } else if (nd[[1]] %in% avoid) {
      search.name[i] <- nms[i]
    } else {
      search.name[i] <- nms[i]
      lng <- .replace(i, 'classification_path')
      if(!is.null(parent)) {
        parent_test <- grepl(parent, lng)
      }
      if(parent_test) {
        name.string[i] <- .replace(i, 'name_string')
        canonical.form[i] <- .replace(i, 'canonical_form')
        lineage[i] <- lng
        lineage.ids[i] <- .replace(i, 'classification_path_ids')
        rank[i] <- .replace(i, 'classification_path_ranks')
        taxid[i] <- .replace(i, 'taxon_id')
        match.type[i] <- .replace(i, 'match_type')
        prescore[i] <- .replace(i, 'prescore')
        score[i] <- nd$results[[1]]$score
      }
    }
  }
  res <- data.frame (search.name=search.name,
                     name.string=name.string,
                     canonical.form=canonical.form,
                     lineage=lineage, lineage.ids=lineage.ids,
                     rank=rank, taxid=taxid,
                     match.type=match.type, prescore=prescore,
                     score=score, stringsAsFactors=FALSE)
  failed <- which (is.na (res$name.string))
  if (genus & length (failed) > 0) {
    # if genus, search just genus names
    genus.nms <- sub ('\\s+.*', '', res$search.name[failed])
    genus.res <- taxaResolve(genus.nms, batch, datasource, genus=FALSE,
                             parent=parent, cache=cache)
    # replace in original results, all slots except search.name
    res[failed,-1] <- genus.res[ ,-1]
  }
  return (res)
}

.safeFromJSON <- function (url, max_trys=5, power=2) {
  # Safe wrapper for fromJSON
  trys <- 0
  waittime <- 2
  while (trys < max_trys) {
    json_obj <- try (RJSONIO::fromJSON(url), silent = TRUE)
    if(class(json_obj) == 'try-error') {
      cat('---- Connection failed: trying again in [', waittime,
          's]----\n', sep='')
      trys <- trys + 1
      Sys.sleep(waittime)
      waittime <- waittime*power
    } else {
      return (json_obj)
    }
  }
  warning("Failed to connect, server may be down.")
  list('data' = NULL)
}

.findClade <- function(lineages) {
  # for a list of lineages, find the clade shared by all
  subj <- lineages[[1]]
  for(i in 2:length(lineages)) {
    query <- lineages[[i]]
    subj <- subj[subj %in% query]
  }
  subj[length(subj)]
}