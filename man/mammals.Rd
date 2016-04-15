\name{mammals}
\alias{mammals}
\docType{data}
\title{Phylogenetic tree of Mammalia}
\description{
  Bininda-Emonds et al. (2007)'s Mammalian supertree augmented with taxonomic information. Here used for testing and demonstrating treeman functions.
}
\usage{
   data(mammals)
}
\format{
  \code{mammals} is a \code{TreeMan} object
}
\source{
   Bininda-Emonds et al. (2007) The Delayed Rise of Present-Day Mammals. Nature, 446(7135): 507-512
}
\examples{
data(mammals) # load object
print(mammals)
print(mammals[['Homo_sapiens']])
}
\keyword{datasets}
