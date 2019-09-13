rm(list=ls())
library(coda)
library(reshape2)
library(ggplot2)
library(SBMSplitMerge)
library(clue)

set.seed(1)

data(Macaque)

ggsave("edges.pdf", image(Macaque))

DMA <- list(edges = edges_bern(), params = param_beta(1,1,1,1), blocks = dma(1, 6))

n_steps <- 15000
sigma   <- 0.1

DMA$output <- sampler(Macaque, DMA, n_steps, "rj", sigma)

burn <- 1000:15000
DMA$plots  <- eval_plots(DMA$output, burn=burn)

ggsave("rj_num_blocks_trace.pdf", DMA$plots$num_blocks_trace)
ggsave("rj_blocks_trace.pdf", DMA$plots$blocks_trace)
ggsave("rj_blocks_trace_sorted.pdf", DMA$plots$blocks_trace_sorted)
ggsave("rj_post_pairs.pdf", DMA$plots$post_pairs)
ggsave("rj_post_pairs_sorted.pdf", DMA$plots$post_pairs_sorted)
ggsave("rj_edges_sorted.pdf", image(edges(Macaque$E[DMA$plots$sortind,DMA$plots$sortind])))

thetas <- t(DMA$output$postt[1,1:5,burn])

## table for latex
DMA$theta_table <- signif(cbind(t(apply(thetas, 2, quantile, probs=c(0.05, 0.5, 0.95))), effectiveSize(thetas)), 2)
colnames(DMA$theta_table) <- c("5%", "50%", "95%", "Effective Size")
rownames(DMA$theta_table) <- paste("theta_", 1:5 -1)

DMA$theta_table

## acceptance
colMeans(diff(thetas) == 0.0)

df <- data.frame(cbind(burn, thetas))
names(df) <- c("index", paste0("theta", 1:ncol(thetas)-1))
mf <- melt(df, id="index")
ggsave("rj_theta.pdf", ggplot(mf, aes(x=index, y=value, col=variable)) + geom_line() + ylim(c(0,1)))

save(DMA, file="rj_output.RData")
