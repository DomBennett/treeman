# take part of the mammalian tree and pin missing tips to it using taxonize

library(treeman)

mammals <- readTree(file="/Users/djb208/Coding/Project-EDBMM/data/raw_trees/literature/bininda.tre", update = FALSE)
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
  nid <- mammals@ndlst[[i]][['id']]
  done <- FALSE
  kids <- getNdKids(mammals, nid)
  if(length(kids) > 0) {
    trms <- gsub("_", " ", kids)
    genus_names <- sapply(strsplit(trms, ' '), function(x) x[[1]])
    if(all(genus_names == genus_names[1])) {
      mammals@ndlst[[i]][['txnym']] <- genus_names[1]
      done <- TRUE
    } else {
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
    }
    if(!done) {
      mammals@ndlst[[i]][['txnym']] <- 'Unknown'
    }
  } else {
    mammals@ndlst[[i]][['txnym']] <- strsplit(mammals@ndlst[[i]][['id']], "_")[[1]][1]
  }
}

mammals@ndlst[['Homo_sapiens']]
mammals@ndlst[['n816']]
mammals@ndlst[['n1044']]

mammalian_orders <- c('Macroscelidea', 'Afrosoricida', 'Tubulidentata', 'Hyracoidea',
                      'Proboscidea', 'Sirenia', 'Pilosa', 'Cingulata', 'Scandentia',
                      'Dermoptera', 'Primates', 'Lagomorpha', 'Rodentia', 'Eulipotyphla',
                      'Cetacea', 'Artiodactyla', 'Chiroptera', 'Perissodactyla', 'Carnivora',
                      'Pholidota')

save(mammals, file="data/mammals.rda", compress="xz")
