library(ggplot2)
library(SBMSplitMerge)
library(xtable)
source("../eval_theta.R")

## big text on plots:
thm <- ggplot2::theme_gray(base_size = 20)

## data
data(Enron)
le <- Enron$Edges$E
le <- edges(ifelse(le > 0, log(le), 0))
ggsave("edges.pdf", image(Enron$Edges) + thm)
ggsave("log_edges.pdf", image(le) + thm)

## params
n_steps <- 10000
n_burns <- 1000
burnin <- (n_burns+1):n_steps
sigma <- 0.5

## RJ NegBin
splitmerge_nbin <- list()
splitmerge_nbin$model  <- sbmmod(dma(1, 10), param_nbin(1,1,1,1,1,1,1,1), edges_nbin())
splitmerge_nbin$output <- sampler(Enron$Edges, splitmerge_nbin$model, n_steps, "rj", sigma, n_steps/10)
splitmerge_nbin$plots <- eval_plots(splitmerge_nbin$output)

range <- 1:10
tmp <- eval_theta(splitmerge_nbin$output$postt, range=range, burnin)
colnames(tmp$theta_table) <- c(matrix(outer(c("r_", "p_"), c("5%", "50%", "95%", "ESS"), paste0),byrow=TRUE))
splitmerge_nbin$theta_table <- tmp$theta_table
splitmerge_nbin$theta_r_plots <- tmp$plots[[1]] +
    scale_color_discrete(name="r", labels=range-1)
splitmerge_nbin$theta_p_plots <- tmp$plots[[2]] +
    scale_color_discrete(name="p", labels=range-1) +
    ylim(c(0,1))
save(splitmerge_nbin, n_steps, sigma, file="rj_nbin.RData")

nbin_le <- log(Enron$Edges$E[splitmerge_nbin$plots$sortind,splitmerge_nbin$plots$sortind])
nbin_le <- ifelse(nbin_le > 0, log(nbin_le), 0)

ggsave("rj_nbin_theta_r_trace.pdf"       , splitmerge_nbin$theta_r_plots + thm)
ggsave("rj_nbin_theta_p_trace.pdf"       , splitmerge_nbin$theta_p_plots + thm)
ggsave("rj_nbin_blocks_trace.pdf"        , splitmerge_nbin$plots$blocks_trace + thm)
ggsave("rj_nbin_blocks_trace_sorted.pdf" , splitmerge_nbin$plots$blocks_trace_sorted + thm)
ggsave("rj_nbin_num_blocks_trace.pdf"    , splitmerge_nbin$plots$num_blocks_trace + thm)
ggsave("rj_nbin_post_pairs.pdf"          , splitmerge_nbin$plots$post_pairs + thm)
ggsave("rj_nbin_post_pairs_sorted.pdf"   , splitmerge_nbin$plots$post_pairs_sorted + thm)
ggsave("rj_nbin_log_edges_sorted.pdf"    , image(edges(nbin_le)) + thm)

## theta table for paper
write.table(signif(splitmerge_nbin$theta_table, 2), "rj_nbin_thetas.tab")


## RJ Pois
splitmerge_pois <- list()
splitmerge_pois$model  <- sbmmod(dma(1, 10), param_gamma(1,1,1,1), edges_pois())
splitmerge_pois$output <- sampler(Enron$Edges, splitmerge_pois$model, n_steps, "rj", sigma, n_steps/10)
splitmerge_pois$plots <- eval_plots(splitmerge_pois$output)

range <- 1:10
tmp <- eval_theta(splitmerge_pois$output$postt, range=range, burnin)
splitmerge_pois$theta_table <- tmp$theta_table
splitmerge_pois$theta_plots <- tmp$plots[[1]] +
    scale_color_discrete(name="lambda", labels=range-1)
save(splitmerge_pois, n_steps, sigma, file="rj_pois.RData")

pois_le <- log(Enron$Edges$E[splitmerge_pois$plots$sortind,splitmerge_pois$plots$sortind])
pois_le <- ifelse(pois_le > 0, log(pois_le), 0)

ggsave("rj_pois_theta_trace.pdf"         , splitmerge_pois$theta_plots + thm)
ggsave("rj_pois_blocks_trace.pdf"        , splitmerge_pois$plots$blocks_trace + thm)
ggsave("rj_pois_blocks_trace_sorted.pdf" , splitmerge_pois$plots$blocks_trace_sorted + thm)
ggsave("rj_pois_num_blocks_trace.pdf"    , splitmerge_pois$plots$num_blocks_trace + thm)
ggsave("rj_pois_post_pairs.pdf"          , splitmerge_pois$plots$post_pairs + thm)
ggsave("rj_pois_post_pairs_sorted.pdf"   , splitmerge_pois$plots$post_pairs_sorted + thm)
ggsave("rj_pois_log_edges_sorted.pdf"    , image(edges(pois_le)) + thm)

## theta table for paper
write.table(signif(splitmerge_pois$theta_table, 2), "rj_pois_thetas.tab")


xtable(splitmerge_nbin$theta_table[,1:3], digits=3)
xtable(splitmerge_nbin$theta_table[,5:7], digits=3)
xtable(splitmerge_pois$theta_table[,1:3], digits=3)
