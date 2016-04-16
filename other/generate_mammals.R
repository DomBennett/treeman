# take part of the mammalian tree and pin missing tips to it using taxonize

library(treeman)

mammals <- readTree(file="other/bininda.tre")
print(mammals)

load("other/mammalia_resolvedlist.Rd")


.findClade <- function (lineages) {
  # for a list of lineages, find the clade shared by all
  subj <- lineages[[1]]
  for (i in 2:length (lineages)) {
    query <- lineages[[i]]
    subj <- subj[subj %in% query]
  }
  subj[length (subj)]
}

for(i in 1:length(mammals@ndlst)) {
  print(mammals@ndlst[[i]][['id']])
  done <- FALSE
  if(!is.null(mammals@ndlst[[i]][['kids']])) {
    trms <- mammals@ndlst[[i]][['kids']]
    trms <- sub("_", " ", trms)
    indxs <- resolve.list$resolved$search.name %in% trms
    if(sum(indxs) > 1) {
      lngs <- resolve.list[['lineages']][indxs]
      bool <- unlist(lapply(lngs, function(l) 'Mammalia' %in% l))
      lngs <- lngs[bool]
      if(length(lngs) > 1) {
        mammals@ndlst[[i]][['txnym']] <- .findClade(lngs)
        done <- TRUE
      }
      if(length(lngs) == 1) {
        lng <- resolve.list[['lineages']][indxs][[1]]
        mammals@ndlst[[i]][['txnym']] <- lng[length(lng)-1]
        done <- TRUE
      }
    }
    if(!done) {
      mammals@ndlst[[i]][['txnym']] <- 'Unknown'
    }
  } else {
    mammals@ndlst[[i]][['txnym']] <- strsplit(mammals@ndlst[[i]][['id']], "_")[[1]]
  }
}

mammals@ndlst[['Homo_sapiens']]
mammals[['n816']]
mammals[['n815']]
mammals[['n814']]
mammals[['n813']]
mammals[['n754']]
mammals[['n753']]
mammals[['n752']]
mammals[['n751']]
mammals[['n750']]
mammals[['n749']]
mammals[['n8']]
mammals[['n7']]
mammals[['n6']]
mammals[['n5']]
mammals[['n4']]
mammals[['n1']]

getTxnym <- function(tree, txnym) {
  # get nd id(s) for a taxonym
  .get <- function(id, txnym, ...) {
    if(rfrnc_txnym %in% txnym) {
      ids <<- c(id, ids)
    }
  }
  rfrnc_txnym <- txnym
  ids <- NULL
  m_ply(tree@ndlst, .fun=.get)
  ids
}

getTxnyms <- function(tree, txnyms) {
  # get nd id(s) for taxonyms
  .get <- function(id, txnym, ...) {
    for(t in txnyms) {
      if(t %in% txnym) {
        res[[t]] <<- c(res[[t]], id)
      }
    }
  }
  res <- list()
  m_ply(tree@ndlst, .fun=.get)
  res
}

mammalian_orders <- c('Macroscelidea', 'Afrosoricida', 'Tubulidentata', 'Hyracoidea',
                      'Proboscidea', 'Sirenia', 'Pilosa', 'Cingulata', 'Scandentia',
                      'Dermoptera', 'Primates', 'Lagomorpha', 'Rodentia', 'Eulipotyphla',
                      'Cetacea', 'Artiodactyla', 'Chiroptera', 'Perissodactyla', 'Carnivora',
                      'Pholidota')

save(mammals, file="data/mammals.rda", compress="xz")
