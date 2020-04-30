library(ggplot2)
library(SBMSplitMerge)
library(coda)
source("../eval_theta.R")

## big text for plots:
thm <- theme_gray(base_size = 20)

## load and plot data
data(Macaque)
ggsave("edges.pdf", image(Macaque) + thm  + xlab("Node") + ylab("Node"))

n_steps <- 10000
n_burn  <- 1000
burnin  <- (n_burn+1):n_steps
sigma   <- 0.1

## RJ
set.seed(1)
splitmerge <- list()
splitmerge$model <- sbmmod(dma(1, 6), param_beta(1,1,1,1), edges_bern())
splitmerge$output <- sampler(Macaque, splitmerge$model, n_steps, "rj", sigma, n_steps/10)
splitmerge$plots  <- eval_plots(splitmerge$output, burnin)
## evaluate thetas - customizing the plot since we have a probability parameter:
range <- 1:5
tmp <- eval_theta(splitmerge$output$postt, range=range, burnin)
splitmerge$theta_table <- tmp$theta_table
splitmerge$theta_plots <- tmp$plots[[1]] + ylim(c(0,1)) +
    scale_color_discrete(name="p", labels=range-1)

write.table(signif(splitmerge$theta_table, 2), "rj_thetas.tab")

save(splitmerge, file="rj_output.RData")
load("rj_output.RData")

ggsave("rj_theta_trace.pdf"         , splitmerge$theta_plots + thm)
ggsave("rj_blocks_trace.pdf"        , splitmerge$plots$blocks_trace + thm)
ggsave("rj_blocks_trace_sorted.pdf" , splitmerge$plots$blocks_trace_sorted + thm)
ggsave("rj_num_blocks_trace.pdf"    , splitmerge$plots$num_blocks_trace + thm)
ggsave("rj_post_pairs.pdf"          , splitmerge$plots$post_pairs + thm)
ggsave("rj_post_pairs_sorted.pdf"   , splitmerge$plots$post_pairs_sorted + thm)
ggsave("rj_edges_sorted.pdf"        , image(edges(Macaque$E[splitmerge$plots$sortind,splitmerge$plots$sortind])) + thm)
##

## ## DP
## set.seed(1)
## dp <- list()
## dp$model  <- sbmmod(crp(5), param_beta(1,1,1,1), edges_bern())
## dp$output <- sampler(Macaque, dp$model, n_steps, "dp", sigma, n_steps/10)
## dp$plots  <- eval_plots(dp$output, burnin)
## ## evaluate thetas - customizing the plot since we have a probability parameter:
## range <- 1:3
## tmp <- eval_theta(dp$output$postt, range=range, burnin)
## dp$theta_table <- tmp$theta_table
## dp$theta_plots <- tmp$plots[[1]] + ylim(c(0,1)) +
##     scale_color_discrete(name="p", labels=range-1)
## save(dp, file="dp_output.RData")
## load("dp_output.RData")
## ggsave("dp_blocks_trace.pdf"        , dp$plots$blocks_trace + thm)
## ggsave("dp_blocks_trace_sorted.pdf" , dp$plots$blocks_trace_sorted + thm)
## ggsave("dp_num_blocks_trace.pdf"    , dp$plots$num_blocks_trace + thm)
## ggsave("dp_post_pairs.pdf"          , dp$plots$post_pairs + thm)
## ggsave("dp_post_pairs_sorted.pdf"   , dp$plots$post_pairs_sorted + thm)
## ggsave("dp_edges_sorted.pdf"        , image(edges(Macaque$E[dp$plots$sortind,dp$plots$sortind])) + thm)
## ##
