library(ggplot2)
library(SBMSplitMerge)

thm <- ggplot2::theme_gray(base_size = 20)
data(Enron)

le <- Enron$Edges$E
le <- edges(ifelse(le > 0, log(le), 0))

ggsave("edges.pdf", image(Enron$Edges) + thm)
ggsave("log_edges.pdf", image(le) + thm)

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
    save(NegBinDMA, n_steps, sigma, file="rj_nbin.RData")
}

{
    PoisDMA <- list(
        blocks = dma(1, 10),
        edges = edges_pois(),
        params = param_gamma(1, 1, 1, 1)
    )
    PoisDMA$rj_output <- sampler(Enron$Edges, PoisDMA, n_steps, "rj", sigma)
    PoisDMA$plots <- eval_plots(PoisDMA$rj_output)
    save(PoisDMA, n_steps, sigma, file="rj_pois.RData")
}

load("rj_nbin.RData")
{
    ## convert edges to log scale and plot sorted by the assignment computed in eval_plots
    nbin_le <- log(Enron$Edges$E[NegBinDMA$plots$sortind,NegBinDMA$plots$sortind])
    nbin_le <- ifelse(nbin_le > 0, log(nbin_le), 0)

    ggsave("rj_nbin_blocks_trace.pdf"        , NegBinDMA$plots$blocks_trace + thm)
    ggsave("rj_nbin_blocks_trace_sorted.pdf" , NegBinDMA$plots$blocks_trace_sorted + thm)
    ggsave("rj_nbin_num_blocks_trace.pdf"    , NegBinDMA$plots$num_blocks_trace + thm)
    ggsave("rj_nbin_post_pairs.pdf"          , NegBinDMA$plots$post_pairs + thm)
    ggsave("rj_nbin_post_pairs_sorted.pdf"   , NegBinDMA$plots$post_pairs_sorted + thm)
    ggsave("rj_nbin_log_edges_sorted.pdf"    , image(edges(nbin_le)) + thm)

    burn <- (NegBinDMA$rj_output$nsteps/2):NegBinDMA$rj_output$nsteps
    inds <- 1:8
    thetas1 <- t(NegBinDMA$rj_output$postt[1,inds,burn])
    thetas2 <- t(NegBinDMA$rj_output$postt[2,inds,burn])

    ## table of metrics
    signif(cbind(t(apply(thetas1, 2, quantile, c(0.05, 0.5, 0.95))), coda::effectiveSize(thetas1)), 2)
    signif(cbind(t(apply(thetas2, 2, quantile, c(0.05, 0.5, 0.95))), coda::effectiveSize(thetas2)), 2)

    df <- data.frame(cbind(burn, thetas1))
    names(df) <- c("index", paste0("theta1_", 1:ncol(thetas1)-1))
    mf <- reshape2::melt(df, id="index")
    ggsave("rj_nbin_theta1.pdf", ggplot(mf, aes(x=index, y=value, col=variable)) + geom_line()+ thm)

    df <- data.frame(cbind(burn, thetas2))
    names(df) <- c("index", paste0("theta2_", 1:ncol(thetas2)-1))
    mf <- reshape2::melt(df, id="index")
    ggsave("rj_nbin_theta2.pdf", ggplot(mf, aes(x=index, y=value, col=variable)) + geom_line()+ thm)
}


load("rj_pois.RData")
{
    ggsave("rj_pois_blocks_trace.pdf"        , PoisDMA$plots$blocks_trace + thm)
    ggsave("rj_pois_blocks_trace_sorted.pdf" , PoisDMA$plots$blocks_trace_sorted + thm)
    ggsave("rj_pois_num_blocks_trace.pdf"    , PoisDMA$plots$num_blocks_trace + thm)
    ggsave("rj_pois_post_pairs.pdf"          , PoisDMA$plots$post_pairs + thm)
    ggsave("rj_pois_post_pairs_sorted.pdf"   , PoisDMA$plots$post_pairs_sorted + thm)

    pois_le <- log(Enron$Edges$E[PoisDMA$plots$sortind,PoisDMA$plots$sortind])
    pois_le <- ifelse(pois_le > 0, log(pois_le), 0)
    ggsave("rj_pois_log_edges_sorted.pdf", image(edges(pois_le)) + thm)

    burn <- (PoisDMA$rj_output$nsteps/2):PoisDMA$rj_output$nsteps
    inds <- 1:8
    thetas <- t(PoisDMA$rj_output$postt[1,inds,burn])
    signif(cbind(t(apply(thetas, 2, quantile, c(0.05, 0.5, 0.95))), coda::effectiveSize(thetas)), 2)

    df <- data.frame(cbind(burn, thetas))
    names(df) <- c("index", paste0("theta1_", 1:ncol(thetas)-1))
    mf <- reshape2::melt(df, id="index")
    ggsave("rj_pois_theta1.pdf", ggplot(mf, aes(x=index, y=value, col=variable)) + geom_line() + thm)
}
