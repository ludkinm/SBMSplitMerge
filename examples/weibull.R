#!/usr/bin/Rscript
library(SBMSplitMerge)

## set up the inference model
model <- list(
    edges = edgemod(
        dg = function(e, p) dweibull(e, p[1,,], p[2,,], log=TRUE)
       ,
        rg = function(p) rweibull(1, p[1], p[2])
    )
   ,
    params = parammod(
        c(1,1,1,1),
        c(1,1,1,1),
        function(n, p)
            cbind(rgamma(n, p[1], p[2]), rgamma(n, p[3], p[4])),
        function(x, p, log)
            dgamma(x[,1], p[1], p[2], log=log) + dgamma(x[,2], p[3], p[4], log=log),
        function(x)
            cbind(log(x[1]), log(x[2])),
        function(x)
            cbind(exp(x[1]), exp(x[2])),
        function(x)
            -log(x[1])-log(x[2])
    )
   ,
    blocks = dma(1, 10)
   ,
    name="weibull_gam_gam"
)

## simulate data
set.seed(1)
true_blocks <- blocks(rep(c(1, 2, 3), c(10, 20, 20)))
true_params <- params(c(1, 0.5), cbind(1, c(3,4,5)))
true_sbm    <- sbm(true_blocks, true_params)
Edges       <- redges(true_params, true_blocks, model$edges)
plot(Edges)

## run sampler
sampler_output <- sampler(Edges, model, 300, "rj", sigma = 0.1)

## evaluate
plots <- eval_plots(sampler_output)
plots$blocks_trace
plots$post_pairs
