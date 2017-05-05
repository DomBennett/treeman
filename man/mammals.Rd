\name{mammals}
\alias{mammals}
\docType{data}
\title{Phylogenetic tree of Mammalia}
\description{
  Bininda-Emonds et al. (2007)'s Mammalian supertree augmented with taxonomic information. Here used for testing and demonstrating treeman functions. See R script to see how the tree was generate: \url{https://github.com/DomBennett/treeman/blob/master/other/generate_tree_data.R}.
}
\usage{
   data(mammals)
}
\format{
  \code{mammals} is a \code{TreeMan} object
}
\source{
   Bininda-Emonds et al. 2007. The Delayed Rise of Present-Day Mammals. Nature, 446(7135): 507-512
}
\examples{
data(mammals) # load object
summary(mammals)
summary(mammals[['Homo_sapiens']])
}
\keyword{datasets}
