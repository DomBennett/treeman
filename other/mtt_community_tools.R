# Tools for 3_turnover.R
# Taken from MoreTreeTools: https://github.com/DomBennett/MoreTreeTools
# For full functionality of functions install MoreTreeTools.

commplot <- function(cmatrix, tree, groups = rep(1, nrow(cmatrix)),
                     no.margin = TRUE, ...){
  # ultrametricize tree FIRST
  tree <- ape::compute.brlen(tree, method="Grafen")
  # plot phylogeny, allow space for points
  tree$edge.length <- tree$edge.length/getSize (tree, 'rtt')
  # make all trees the same length before plotting i.e. all branches from terminal
  # node to tip equal 1
  # for some weird reason the rules of plotting are dyanmic!
  if (nrow(cmatrix) < 20) {
    variable.max <- 1 + nrow(cmatrix)/20
    spacing.start <- 0.55
    spacing.i <- 0.05
  } else {
    variable.max <- nrow(cmatrix)/10
    spacing.i <- 0.1 - 1/nrow(cmatrix)
    spacing.start <- 0.5 + spacing.i
  }
  plot(tree, no.margin = no.margin, show.tip.label = FALSE,
       x.lim = c(0, variable.max), ...)
  # generate alphas based on abundances
  n <- length(unique(groups))
  hs <- seq.int(0, 1 + max(1, n - 1)/n, length.out = n)%%1
  alphas <- cmatrix/max(cmatrix)
  # loop init
  ntips <- length(tree$tip.label)
  spacing <- spacing.start
  group <- groups[1]
  j <- 1
  # loop through sites and plot points for present species
  for(i in 1:nrow(cmatrix)){
    j <- ifelse(group == groups[i], j, j + 1)
    pull <- as.logical(cmatrix[i,])
    taxa <- tree$tip.label[pull]
    abunds <- alphas[i, pull]
    ape::tiplabels(tip = match(taxa, tree$tip.label),
              pch = 19, adj = spacing, col = hsv(rep(hs[j], ntips), 1, 1, abunds))
    spacing <- spacing + spacing.i
    group <- groups[i]
  }
}

genCommData <- function(tree, focal, psi = 1, mean.incid, mean.abun = FALSE,
                        nsites = 1) {
  ##TODO: write test
  invertVector <- function(dists) {
    # for reversing the probs for overdispersion
    u <- sort(unique(dists), TRUE)
    s <- sort(u)
    probs <- rep(NA, length(dists))
    for (i in 1:length(u)) {
      probs[u[i] == dists] <- s[i]
    }
    return (probs)
  }
  genAbuns <- function(row) {
    # for generating abundances for each row
    out.row <- rep(0, ntips)
    temp.probs <- probs
    temp.probs[row < 1] <- 0
    abundance <- abs(ceiling(rnorm(1, mean = mean.abun)))
    if (abundance == 0) {
      return (out.row)
    } else {
      abuns <- sample(1:ntips, abundance, prob = temp.probs, replace = TRUE)
      abuns <- table(abuns)
      out.row[as.numeric(names(abuns))] <- abuns
      return (out.row)
    }
  }
  ntips <- length(tree$tip.label)
  output <- matrix(rep(NA, ntips * nsites),
                   ncol = ntips, nrow = nsites)
  colnames(output) <- tree$tip.label
  pd.dists <- ape::cophenetic.phylo(tree)
  focal.dists <- pd.dists[ , focal] + 1 # avoid 0
  if (psi > 0) focal.dists <- invertVector(focal.dists)
  probs <- focal.dists^abs(psi)
  probs <- probs/sum(probs)
  for (i in 1:nsites) {
    incidence <- abs(ceiling(rnorm(1, mean = mean.incid))) # avoid negative numbers
    output[i, ] <- ifelse(1:ntips %in% sample(1:ntips, incidence,
                                              prob = probs), 1, 0)
  }
  if (mean.abun != FALSE) {
    for (i in 1:nsites) {
      output[i, ] <- genAbuns(output[i, ])
    }
  }
  return (output)
}

getSize <- function (tree, type = c ('ntips', 'pd', 'rtt')) {
  type = match.arg (type)
  if (type == 'ntips') {
    return (length (tree$tip.label))
  } else if (type == 'pd') {
    return (sum (tree$edge.length))
  } else {
    age <- max (diag (ape::vcv.phylo (tree)))
    if (!is.null (tree$root.edge)) {
      age <- age + tree$root.edge
    }
    return (age)
  }
}