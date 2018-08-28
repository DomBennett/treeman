
# TODO: 
# -- rework how upgrading works
# -- reGravitise and reBalance
# -- create node object methods
# -- save and load TreeMen


library(treeman)
txt<-"(Mimus_gilvus:42.797302,((((Coereba_flaveola:8.469427,Tiaris_bicolor:8.469427)Node5:8.287735,Sicalis_flaveola:16.757164)Node4:4.542213,(Setophaga_ruticilla:19.005642,(Zonotrichia_capensis:16.089878,(Quiscalus_major:13.013371,(Icterus_nigrogularis:8.389537,Icterus_icterus:8.389537)Node9:4.623834)Node8:3.076508)Node7:2.915765)Node6:2.293731)Node3:8.74406,Passer_domesticus:30.043432)Node2:12.753868)Node1:73.14138;"
phylo <- ape::read.tree(text=txt)
tree <- as(phylo, 'TreeMan')

summary(tree)
