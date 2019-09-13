#!/usr/bin/Rscript
rm(list=ls())

library(ggplot2)
library(SBMSplitMerge)

cmdargs <- commandArgs(TRUE)
seed    <- get_arg(cmdargs, 1, 1)

set.seed(seed)
true_blocks <- blocks(rep(c(1, 2, 3, 4), c(19, 23, 27, 31)))
true_params <- params(c(1, 0.5), cbind(c(1,4:6), 0.5))
true_sbm    <- sbm(true_blocks, true_params)
Edges <- redges(true_params, true_blocks, edges_nbin())

save(Edges, file="nbin/edges.RData")
ggsave(paste0("nbin/edges.pdf"), image(Edges))
ggsave(paste0("nbin/edges_blocks.pdf"), image(Edges, true_sbm))
