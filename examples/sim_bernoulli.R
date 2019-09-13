#!/usr/bin/Rscript
rm(list=ls())

library(ggplot2)
library(SBMSplitMerge)

cmdargs <- commandArgs(TRUE)
seed    <- get_arg(cmdargs, 1, 1)

set.seed(seed)
true_blocks <- blocks(rep(c(1, 2, 3, 4), c(19, 23, 27, 31)))
true_params <- params(0.05, (3+1:4)/10)
true_sbm    <- sbm(true_blocks, true_params)
Edges <- redges(true_params, true_blocks, edges_bern())

save(Edges, file="bern/edges.RData")
ggsave(paste0("bern/edges.pdf"), image(Edges))
ggsave(paste0("bern/edges_blocks.pdf"), image(Edges, true_sbm))
