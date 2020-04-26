library(SBMSplitMerge)

## set up the model
model <- sbmmod(dma(1, 10), param_norm(0,5,1,1,2,5,1,1), edges_norm(), name="norm")

## simulate data
set.seed(1)
true_blocks <- blocks(rep(c(1, 2, 3, 4), c(19, 23, 27, 31)))
true_params <- params(c(0, 1), cbind(c(3,4,5,6), 1))
true_sbm    <- sbm(true_blocks, true_params)
Edges       <- redges(true_sbm, model$edge)
plot(Edges)

## run sampler
sampler_output <- sampler(Edges, model, 100, "rj", 0.5)

## evaluate
plots <- eval_plots(sampler_output)
plots$post_pairs
