# treeman
[![Build Status](https://travis-ci.org/DomBennett/treeman.svg)](https://travis-ci.org/DomBennett/treeman)[![Coverage Status](https://coveralls.io/repos/DomBennett/treeman/badge.svg?branch=master&service=github)](https://coveralls.io/github/DomBennett/treeman?branch=master)

> An R package for manipulating phylogentic trees using an intuitive S4 class structure.

The `treeman` R package provides a `list` based class for encoding phylogenetic trees in R, making manipulating phylogenetic trees easier to code and more efficient to run.

**Design features**

* Lightweight (few dependencies)
* Fast (vectorised or recursive)
* Simple and intuitive

**Quick guide**

```{R}
# install development copy
library(devtools)
install_github('dombennett/treeman')
# working with the TreeMan class
library(treeman)
?TreeMan  # check the documentation
tree <- randTree(10)  # generate a random tree of 10 tips
print(tree)  # check key stats
tree["tips"]  # extract key stats
```

For more details check out the [wiki](https://github.com/DomBennett/treeman/wiki).

**Licence**

GPL-2

**Status**

In development

**Author**

D.J. Bennett (*but I welcome pull requests!*)
