\name{plants}
\alias{plants}
\docType{data}
\title{Phylogenetic tree of European Flora}
\description{
Durka and Michalski (2012)'s large phylogenetic tree of European vascular plants augmented with taxonomic information generated from the NCBI taxonomy. Here used for testing and demonstrating treeman functions. See R script to see how the tree was generate: \url{https://github.com/DomBennett/treeman/blob/master/other/generate_tree_data.R}.
}
\usage{data(plants)}
\format{\code{plants} is a \code{TreeMan} object}
\source{
    Durka, W. and Michalski, S.G. 2012. Daphne: a dated phylogeny of a large European flora for phylogenetically informed ecological analyses. Ecology 93:2297-2297.
}
\examples{
data(plants) # load object
summary(plants)
}
\keyword{datasets}