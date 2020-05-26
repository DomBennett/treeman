\name{treeman-package}
\alias{treeman-package}
\alias{treeman}
\docType{package}
\title{
Phylogenetic Tree Manipulation
}
\description{
Manipulate phylogenetic trees in R simply, intuitively and efficiently with
a list-based tree structure. Supported methods include: adding and removing
tips or nodes, rescaling branches, calculating common phylogenetic metrics plus
more!
}
\details{
\tabular{ll}{
Package: \tab treeman\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-12-02\cr
License: \tab GPL-2\cr
}
}
\author{
D.J. Bennett

Maintainer: D.J. Bennett <dominic.john.bennett@gmail.com>
}
\references{
Bennett, D.J., Sutton, M.D. & Turvey, S.T., 2017. treeman: an R package for efficient and intuitive manipulation of phylogenetic trees. BMC Research Notes, 10(1), p.30. DOI: 10.1186/s13104-016-2340-8
}
\concept{phylogeny}
\concept{evolution}
\concept{macroevolution}
\concept{tree simulation}
\seealso{
\url{https://github.com/DomBennett/treeman/wiki}
}
\examples{
library(treeman)
tree <- randTree(100)
print(tree)
}
