# treeman

An R package for manipulating phylogentic trees using a simple phylogenetic S4 class structure.

## Background

Adding and removing tips can be very slow for large trees, this can make processes
that require a tree to be manipulated (e.g. simulation, collapsing, merging) inhibitively
slow. The idea for \code{treeman} is to provide simple \code{list} representation of a
tree that means adding and removing tips or clades is computationally equivalent to
adding and removing elements to a list. This also allows manipulations to be readily vectorised.

## Goals

* Lightweight (few dependencies)
* Fast (vectorised or recursive)
* Tested
* Simple