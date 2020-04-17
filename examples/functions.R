library(ggplot2)

## run a single chain producing output plots and saving to disk
Run <- function(Edges, model, n_steps, n_burn, algorithm="rj", sigma=0.1){
    fname <- paste0(model$name, "/", algorithm, "_post_", n_steps)
    if(missing(n_burn))
        n_burn <- n_steps/2
    sampler_output <- sampler(Edges, model, n_steps, algorithm, sigma,
                              currsbm=rsbm(Edges$numnodes, model))
    plots <- eval_plots(sampler_output, n_burn1:n_steps)
    post_theta <- sampler_output$postt[,,n_burn:n_steps, drop=FALSE]
    post_theta_summary <- array(dim=c(dim(post_theta)[1], 3, dim(post_theta)[2]))
    for(i in 1:dim(post_theta)[1])
        post_theta_summary[i,,] <- apply(post_theta[i,,], 1, quantile, probs=c(.05,.5,.95), na.rm=TRUE)
    save(algorithm, sampler_output, plots, Edges, model,
         n_steps, sigma, post_theta_summary,
         file=paste0(fname, ".RDa"))
    ggsave(file=paste0(fname , "_num_blocks_trace.pdf")    , plots$num_blocks_trace    + theme_gray(base_size = 20))
    ggsave(file=paste0(fname , "_blocks_trace.pdf")        , plots$blocks_trace        + theme_gray(base_size = 20))
    ggsave(file=paste0(fname , "_blocks_trace_sorted.pdf") , plots$blocks_trace_sorted + theme_gray(base_size = 20))
    ggsave(file=paste0(fname , "_param_trace.pdf") ,           plots$param_trace + theme_gray(base_size = 20))
    ggsave(file=paste0(fname , "_param_trace_sorted.pdf") , plots$param_trace_sorted + theme_gray(base_size = 20))
    ggsave(file=paste0(fname , "_post_pairs.pdf") , plots$post_pairs + theme_gray(base_size = 20))
    ggsave(file=paste0(fname , "_post_pairs_sorted.pdf") , plots$post_pairs_sorted + theme_gray(base_size = 20))
    post_theta_summary
}


## run two chains, one starting with all nodes in block, another with all nodes in a block each
## plot the trace of number of blocks to show they converge to the same thing
PerfectSimulation <- function(Edges, model, n_steps, algorithm="rj", sigma=0.1){
    ## all in one block
    lower <- sampler(Edges, model, n_steps, algorithm, sigma,
                     currsbm=sbm(blocks(rep(1, Edges$numnodes)), rparams(1, model$params)))
    ## all in seperate blocks
    upper <- sampler(Edges, model, n_steps, algorithm, sigma,
                     currsbm=sbm(blocks(1:Edges$numnodes), rparams(Edges$numnodes, model$params)))
    save(lower, upper, Edges, model, n_steps, algorithm, sigma,
         file=paste0(model$name, "/", algorithm, "_perfect_post_", n_steps, ".RDa"))
    ## plots
    perfect <- data.frame(x=1:n_steps, upper=upper$postk, lower=lower$postk)
    save(perfect, file=paste0(model$name, "/", algorithm, "_perfect_",  n_steps, ".RDa"))
    ggsave(paste0(model$name, "/", algorithm, "_perfect_trace_",  n_steps, ".pdf"),
           ggplot(perfect, aes(x, upper, lower)) +
           geom_line(aes(y = lower, colour = "lower")) +
           geom_line(aes(y = upper, colour = "upper")) +
           xlab("Step") + ylab("Number of blocks") +
           labs(color="Chain") +
           scale_y_log10() + theme_gray(base_size = 20),
           width=5, height=5)
    perfect
}

## parallel run n_chains chains each for 2*n_steps steps and compute gelman-rubin diagnostics
GelmanRubin <- function(Edges, model, n_steps, algorithm="rj", n_chains=30, sigma=0.1){
    chains <- parallel::mclapply(1:n_chains, function(n_chains) sampler(Edges, model, 2*n_steps, algorithm, sigma), mc.cores=n_chains)
    x <- matrix(NA, n_chains, n_steps)
    y <- matrix(NA, n_chains, n_steps)
    thetaMeans <- list()
    thetaVars <- list()
    for(j in 1:n_chains){
        for(i in 1:n_steps){
            if(!is.null(chains[[j]])){
                x[j,i] <- mean(chains[[j]]$postt[1,,n_steps+i], na.rm=TRUE)
                y[j,i] <- var(chains[[j]]$postt[1,, n_steps+i], na.rm=TRUE)
            }
            thetaMeans[[j]] <- coda::as.mcmc(x[j,])
            thetaVars[[j]]  <- coda::as.mcmc(y[j,])
        }
    }
    means <- coda::gelman.diag(thetaMeans)
    vars  <- coda::gelman.diag(thetaVars)
    save(means, vars,
         file=paste0(model$name, "/", algorithm, "_gelman_rubin_",  n_steps, ".RDa"))
    list(means=means, vars=vars)
}
