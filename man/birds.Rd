\name{birds}
\alias{birds}
\docType{data}
\title{Phylogenetic tree of Aves}
\description{
  Jetz et al. (2012)'s Avian supertree augmented with taxonomic information generated from the NCBI taxonomy. Here used for testing and demonstrating treeman functions. See R script to see how the tree was generate: \url{https://github.com/DomBennett/treeman/blob/master/other/generate_tree_data.R}.
}
\usage{
   data(birds)
}
\format{
  \code{birds} is a \code{TreeMan} object
}
\source{
  Jetz, W., Thomas, G.H., Joy, J.B., Hartmann, K. and Mooers, A.O. 2012. The global diversity of birds in space and time. Nature, 491: 444-448
}
\examples{
data(birds) # load object
summary(birds)
}
\keyword{datasets}
