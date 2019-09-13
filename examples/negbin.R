#!/usr/bin/Rscript
rm(list=ls())
library(coda)
library(ggplot2)
library(parallel)
library(SBMSplitMerge)

source("functions.r")

cmdargs <- commandArgs(TRUE)

n_steps  <- get_arg(cmdargs, 1, 100)
sigma    <- get_arg(cmdargs, 2, 0.1)
seed     <- get_arg(cmdargs, 3, 1)
n_chains <- get_arg(cmdargs, 4, 1)

model <- list(edges = edges_nbin(), params = param_nbin(1, 1, 1, 1, 1, 1, 1, 1), name="nbin")

load(paste0(model$name, "/edges.RData"))

set.seed(seed)
model$blocks <- dma(1, 10)
runRJ(n_steps, Edges, model, n_chains)

set.seed(seed)
model$blocks <- crp(10)
runDP(n_steps, Edges, model, n_chains)



load("./nbin/rj_post_10000.RDa")
signif(apply(t(rj_output$postt[1,1:5,5000:10000]), 2, quantile, probs=c(.05,.5,.95)), 3)
signif(apply(t(rj_output$postt[2,1:5,5000:10000]), 2, quantile, probs=c(.05,.5,.95)), 3)

load("./nbin/rj_rubin_gelman.RDa")
RJmeans$psrf
RJvars$psrf
