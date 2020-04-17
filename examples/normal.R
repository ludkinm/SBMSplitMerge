#!/usr/bin/Rscript
library(SBMSplitMerge)

## set up the inference model
model <- list(edges = edges_norm(),
              params = param_norm(1, 1, 1, 1, 1, 1, 1, 1),
              blocks = dma(1, 10),
              name="norm")

## simulate data
set.seed(1)
true_blocks <- blocks(rep(c(1, 2, 3, 4), c(19, 23, 27, 31)))
true_params <- params(c(0, 1), cbind(c(3,4,5,6), 1))
true_sbm    <- sbm(true_blocks, true_params)
Edges       <- redges(true_params, true_blocks, model$edges)
plot(Edges)

## run sampler
sampler_output <- sampler(Edges, model, 100, "rj", 0.5)

## evaluate
plots <- eval_plots(sampler_output)
plots$post_pairs
