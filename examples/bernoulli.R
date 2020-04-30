library(SBMSplitMerge)
source("functions.R")

dir.create("bern")

## set up the model
model <- sbmmod(dma(1, 10), param_beta(1, 1, 1, 1), edges_bern(), name="bern")

## simulate data
set.seed(1)
true_blocks <- blocks(rep(c(1, 2, 3, 4), c(19, 23, 27, 31)))
true_params <- params(0.05, (3+1:4)/10)
true_sbm    <- sbm(true_blocks, true_params)
Edges       <- redges(true_sbm, model$edge)
save(Edges, file="bern/edges.RData")
save(true_blocks, true_params, true_sbm, file="bern/truth.RData")
ggsave(paste0("bern/edges.pdf"), image(Edges) + theme_gray(base_size = 20) + xlab("Node") + ylab("Node"))
ggsave(paste0("bern/edges_blocks.pdf"), image(Edges, true_sbm) + theme_gray(base_size = 20) + xlab("Node") + ylab("Node"))

## run a single chain
set.seed(1)
n_steps <- 10000
sigma   <- 0.5
post_theta_summary <- Run(Edges, model, n_steps, sigma=sigma, statusfreq=n_steps/10)

## perfect simulation diagnostic
set.seed(1)
PerfectSimulation(Edges, model, n_steps, sigma=sigma)

## compute gelman-rubin statistics
set.seed(1)
GelmanRubin(Edges, model, n_steps, sigma=sigma)

load("./bern/rj_post_10000.RDa")
tmp <- eval_plots(sampler_output, theta_index = 1:9)
tmp$param_trace[[1]] + geom_hline(yintercept = true_params$theta0) + geom_hline(yintercept = true_params$thetak)
