# treeman

An R package for manipulating phylogentic trees using a simple phylogenetic S4 class structure.

## Background

Adding and removing tips can be very slow for large trees, this can make processes
that require a tree to be manipulated (e.g. simulation, collapsing, merging) inhibitively
slow. The idea for `treeman` is to provide a simple `list` representation of a
tree that means adding and removing tips or clades is computationally equivalent to
adding and removing elements to a list. This also allows manipulations to be readily vectorised.

## Goals

* Lightweight (few dependencies)
* Fast (vectorised or recursive)
* Tested
* Simple

## Quick guide
```{R}
# install
library(devtools)
install_github('dombennett/treeman')
# working with the TreeMan class
library(treeman)
?TreeMan
tree <- randTree(10)
print(tree)
```
```{bash}
Tree (TreeMan Object):
  -- [10] tips
  -- [9] internal nodes
  -- Binary
  -- Age [3.49]
  -- PD [7.66]
  -- Not ultrametric (with extinct tips)
  -- Root node is ["n1"]
```
```{R}
# simple visualization
viz(tree)
```
![tree-viz](https://raw.githubusercontent.com/DomBennett/treeman/master/other/viz-tree.jpeg)

(See [ggtree](https://github.com/GuangchuangYu/ggtree) for more advanced plotting.)

## Status
In development

## Author
D.J. Bennett