# LIBS
library(treeman)
library(testthat)

# DATA
test.lineages <- 
  list (c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyA', 'genusA',
           'speciesA'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyA', 'genusA',
           'speciesB'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyA', 'genusA',
           'speciesC'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyA', 'genusA',
           'speciesD'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyA', 'genusB',
           'speciesE'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyA', 'genusB',
           'speciesF'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyA', 'genusC',
           'speciesG'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyB', 'genusD',
           'speciesH'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyC', 'genusE',
           'speciesI'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderB', 'familyD', 'genusF',
           'speciesJ'))

# RUNNING
context('Testing \'server-methods\'')
test_that ('.safeFromJSON([basic]) works', {
  # simply show that an error is thrown and handled
  expect_that (
    treeman:::.safeFromJSON (
      url='dummyaddress', max_trys=0), throws_error())
})
test_that ('taxaResolve([basic]) works', {
  test_names <- c("Macaca mulatta", "Gorilla gorilla", "Homo sapiens", "Pan paniscus",
                  "Pan troglodytes", "Pongo pygmaeus", "Hylobates agilis",
                  "Hylobates lar", "Hylobates moloch", "Hylobates muelleri", 
                  "Hylobates pileatus", "Hylobates klossii", "Hylobates hoolock",
                  "Hylobates syndactylus", "Hylobates concolor", "Hylobates leucogenys",
                  "Hylobates gabriellae", 'thisisnotaname')
  expected_dimensions <- c (17, 10)
  res <- taxaResolve(test_names)
  res <- res[complete.cases(res), ]
  expect_that(dim(res), equals(expected_dimensions))
})
test_that ('searchTxnyms([basic]) works', {
  tree <- randTree(8)
  new_tids <- c("Gallus_gallus", "Aileuropoda_melanoleucha", "Ailurus_fulgens",
                "Rattus_rattus", "Mus_musculus", "Gorilla_gorilla", "Pan_trogoldytes", "Homo_sapiens")
  tree <- setNdsID(tree, tree['tips'], new_tids)
  tree <- updateTree(tree)
  nd_labels <- searchTxnyms(tree)
  expect_that(sum(is.na(nd_labels)), equals(0))
})
test_that ('.findClade([basic]) works', {
  # class A is shared by all test species
  expect_that (treeman:::.findClade (
    test.lineages), equals ('classA'))
})