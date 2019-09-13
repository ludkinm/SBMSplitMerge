## code to prepare `enron` dataset goes here

library(usethis)
library(igraphdata)
library(SBMSplitMerge)

data(enron)
E <- igraph::as_adjacency_matrix(enron)
Edges <- edges(as.matrix(E))
Group <- igraph::get.vertex.attribute(enron)$Note
Enron <- list(Edges=Edges, Group=Group)
usethis::use_data(Enron)
