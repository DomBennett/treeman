# treeman
[![Build Status](https://travis-ci.org/DomBennett/treeman.svg)](https://travis-ci.org/DomBennett/treeman)[![Coverage Status](https://coveralls.io/repos/DomBennett/treeman/badge.svg?branch=master&service=github)](https://coveralls.io/github/DomBennett/treeman?branch=master)[![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/treeman)](https://cran.r-project.org/web/packages/treeman/index.html)

> An R package for manipulating phylogentic trees using an intuitive S4 class structure.

The `treeman` R package provides a `list` based class for encoding phylogenetic trees in R, making manipulating phylogenetic trees easier to code and more efficient to run.

**Design features**

* Lightweight (few dependencies)
* Fast (vectorised or recursive)
* Simple and intuitive

**Installation**

With CRAN:

```{R}
install.packages('treeman')
```

Installing the development copy via GitHub
```{R}
library(devtools)
install_github('dombennett/treeman')
```

**Quick guide**

```{R}
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

Version 1 released.

**Author**

D.J. Bennett (*but I welcome pull requests!*)
