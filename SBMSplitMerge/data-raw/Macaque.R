## code to prepare `Macaque` dataset goes here

library(usethis)
library(igraphdata)
library(SBMSplitMerge)

data(macaque)
E <- igraph::as_adjacency_matrix(macaque)
Macaque <- edges(as.matrix(E))
usethis::use_data(Macaque)
