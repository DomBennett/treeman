# treeman
[![Build Status](https://travis-ci.org/DomBennett/treeman.svg)](https://travis-ci.org/DomBennett/treeman)[![Coverage Status](https://coveralls.io/repos/DomBennett/treeman/badge.svg?branch=master&service=github)](https://coveralls.io/github/DomBennett/treeman?branch=master)[![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/treeman)](https://CRAN.R-project.org/package=treeman)[![Rdoc](http://www.rdocumentation.org/badges/version/treeman)](http://www.rdocumentation.org/packages/treeman) 

> An R package for manipulating phylogentic trees using an intuitive S4 class structure.

The `treeman` R package provides a `list` based class for encoding phylogenetic trees in R, making manipulating phylogenetic trees easier to code and more efficient to run. `treeman` is aimed to be fast, simple and intuitive.

**Installation**

With CRAN:

```r
install.packages('treeman')
```

Installing the development copy via GitHub
```r
library(devtools)
install_github('dombennett/treeman')
```

**Quick guide**

```r
# working with the TreeMan class
library(treeman)
?TreeMan  # check the documentation
tree <- randTree(10)  # generate a random tree of 10 tips
summary(tree)  # check key stats
tree["tips"]  # extract key stats
```

For more details check out the [wiki](https://github.com/DomBennett/treeman/wiki).

**Licence**

GPL-2

**Status**

Version 1.1 released.

**Reference**

Bennett, D.J., Sutton, M.D. & Turvey, S.T., 2017. treeman: an R package for efficient and intuitive manipulation of phylogenetic trees. *BMC Research Notes*, 10(1), p.30. [Available online](https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-2340-8)

**Author**

D.J. Bennett (*but I welcome pull requests!*)
