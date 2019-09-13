rm(list=ls())
library(coda)
library(ggplot2)
library(reshape2)
library(SBMSplitMerge)

data(Enron)

le <- Enron$Edges$E
le <- ifelse(le > 0, log(le), 0)

ggsave("edges.pdf", image(Enron$Edges))
ggsave("log_edges.pdf", image(le))

n_steps <- 10000
sigma <- 0.5

{
    NegBinDMA <- list(
        blocks = dma(1, 10),
        edges = edges_nbin(),
        params = param_nbin(1, 1, 1, 1,
                            1, 1, 1, 1)
    )
    NegBinDMA$rj_output <- sampler(Enron$Edges, NegBinDMA, n_steps, "rj", sigma)
    NegBinDMA$plots <- eval_plots(NegBinDMA$rj_output)
}


{
    PoisDMA <- list(
        blocks = dma(1, 10),
        edges = edges_pois(),
        params = param_gamma(1, 1, 1, 1)
    )
    PoisDMA$rj_output <- sampler(Enron$Edges, PoisDMA, n_steps, "rj", sigma)
    PoisDMA$plots <- eval_plots(PoisDMA$rj_output)
}

{
    ggsave("rj_nbin_num_blocks_trace.pdf", NegBinDMA$plots$num_blocks_trace)
    ggsave("rj_nbin_blocks_trace.pdf", NegBinDMA$plots$blocks_trace)
    ggsave("rj_nbin_blocks_trace_sorted.pdf", NegBinDMA$plots$blocks_trace_sorted)
    ggsave("rj_nbin_post_pairs.pdf", NegBinDMA$plots$post_pairs)
    ggsave("rj_nbin_post_pairs_sorted.pdf", NegBinDMA$plots$post_pairs_sorted)

    nbin_le <- log(Enron$Edges$E[NegBinDMA$plots$sortind,NegBinDMA$plots$sortind])
    nbin_le <- ifelse(nbin_le > 0, log(nbin_le), 0)
    image(edges(nbin_le))

    ggsave("rj_nbin_log_edges_sorted.pdf", image(edges(nbin_le)))

    burn <- (n_steps/2):n_steps
    inds <- 1:8
    thetas1 <- t(NegBinDMA$rj_output$postt[1,inds,burn])
    thetas2 <- t(NegBinDMA$rj_output$postt[2,inds,burn])

    ## tables
    signif(cbind(t(apply(thetas1, 2, quantile, c(0.05, 0.5, 0.95))), effectiveSize(thetas1)), 2)
    signif(cbind(t(apply(thetas2, 2, quantile, c(0.05, 0.5, 0.95))), effectiveSize(thetas2)), 2)

    df <- data.frame(cbind(burn, thetas1))
    names(df) <- c("index", paste0("theta1_", 1:ncol(thetas1)-1))
    mf <- melt(df, id="index")
    ggsave("rj_nbin_theta1.pdf", ggplot(mf, aes(x=index, y=value, col=variable)) + geom_line())

    df <- data.frame(cbind(burn, thetas2))
    names(df) <- c("index", paste0("theta2_", 1:ncol(thetas2)-1))
    mf <- melt(df, id="index")
    ggsave("rj_nbin_theta2.pdf", ggplot(mf, aes(x=index, y=value, col=variable)) + geom_line())

    save(NegBinDMA, n_steps, sigma, file="rj_nbin.RData")
}

{
    ggsave("rj_pois_num_blocks_trace.pdf", PoisDMA$plots$num_blocks_trace)
    ggsave("rj_pois_blocks_trace.pdf", PoisDMA$plots$blocks_trace)
    ggsave("rj_pois_blocks_trace_sorted.pdf", PoisDMA$plots$blocks_trace_sorted)
    ggsave("rj_pois_post_pairs.pdf", PoisDMA$plots$post_pairs)
    ggsave("rj_pois_post_pairs_sorted.pdf", PoisDMA$plots$post_pairs_sorted)

    pois_le <- log(Enron$Edges$E[PoisDMA$plots$sortind,PoisDMA$plots$sortind])
    pois_le <- ifelse(pois_le > 0, log(pois_le), 0)
    image(edges(pois_le))

    ggsave("rj_pois_log_edges_sorted.pdf", image(edges(pois_le)))

    burn <- (n_steps/2):n_steps
    inds <- 1:8
    thetas <- t(PoisDMA$rj_output$postt[1,inds,burn])
    signif(cbind(t(apply(thetas, 2, quantile, c(0.05, 0.5, 0.95))), effectiveSize(thetas)), 2)

    df <- data.frame(cbind(burn, thetas))
    names(df) <- c("index", paste0("theta1_", 1:ncol(thetas)-1))
    mf <- melt(df, id="index")
    ggsave("rj_pois_theta1.pdf", ggplot(mf, aes(x=index, y=value, col=variable)) + geom_line())

    save(PoisDMA, n_steps, sigma, file="rj_pois.RData")
}



load("rj_nbin.RData")

NegBinDMA$plots$sortind
NegBinDMA$rj_output$postz[72,]
