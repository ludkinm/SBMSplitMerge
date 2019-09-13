## For RJMCMC split-merge sampler
runRJ <- function(n_steps, Edges, model, n_chains=1){
    if(n_chains == 1){
        init_sbm <- rsbm(Edges$numnodes, model)
        rj_output <- sampler(Edges, model, n_steps, "rj", sigma, currsbm=init_sbm)
        save(rj_output, model, Edges, sigma, n_steps, file=paste0(model$name, "/rj_post_", n_steps, ".RDa"))
        rj_plots <- eval_plots(rj_output, 1:n_steps)
        for(i in 1:length(rj_plots))
            ggsave(file=paste0(model$name, "/rj_", names(rj_plots)[i], ".pdf"), rj_plots[[i]])
    } else if(n_chains == 0){
        ## perfect simulation - start at all in one or all in different block configs
        ## all in one block
        init_sbm <- sbm(blocks(rep(1, Edges$numnodes)), rparams(1, model$params))
        lower    <- sampler(Edges, model, n_steps, "rj", sigma, currsbm=init_sbm)
        ## all in seperate blocks
        init_sbm <- sbm(blocks(1:Edges$numnodes), rparams(Edges$numnodes, model$params))
        upper    <- sampler(Edges, model, n_steps, "rj", sigma, currsbm=init_sbm)
        ## plots
        perfectRJ <- data.frame(x=1:n_steps, upper=upper$postk, lower=lower$postk)
        save(perfectRJ, file=paste0(model$name, "/rj_perfect.RDa"))
        ggsave(paste0(model$name, "/rj_perfect_trace.pdf"),
               ggplot(perfectRJ, aes(x, upper, lower)) +
               geom_line(aes(y = lower, colour = "lower")) +
               geom_line(aes(y = upper, colour = "upper")) +
               xlab("Iteration") + ylab("Number of blocks") +
               labs(color="Chain") +
               scale_y_log10(),
               width=5, height=5)
    } else{
        RJ <- mclapply(1:n_chains, function(n_chains) sampler(Edges, model, 2*n_steps, "rj", sigma), mc.cores=n_chains)
        x <- matrix(NA, n_chains, n_steps)
        y <- matrix(NA, n_chains, n_steps)
        thetaMeans <- list()
        thetaVars <- list()
        for(j in 1:n_chains){
            for(i in 1:n_steps){
                if(!is.null(RJ[[j]])){
                    x[j,i] <- mean(RJ[[j]]$postt[1,,n_steps+i], na.rm=TRUE)
                    y[j,i] <- var(RJ[[j]]$postt[1,, n_steps+i], na.rm=TRUE)
                }
                thetaMeans[[j]] <- as.mcmc(x[j,])
                thetaVars[[j]] <- as.mcmc(y[j,])
            }
        }
        RJmeans <- gelman.diag(thetaMeans)
        RJvars  <- gelman.diag(thetaVars)
        save(RJmeans, RJvars, file=paste0(model$name, "/rj_rubin_gelman.RDa"))
    }
}

## For RJMCMC split-merge sampler
runDP <- function(n_steps, Edges, model, n_chains=1){
    if(n_chains == 1){
        init_sbm <- rsbm(Edges$numnodes, model)
        dp_output <- sampler(Edges, model, n_steps, "dp", sigma, currsbm=init_sbm)
        save(dp_output, model, Edges, sigma, n_steps, file=paste0(model$name, "/dp_post_", n_steps, ".RDa"))
        dp_plots <- eval_plots(dp_output, 1:n_steps)
        for(i in 1:length(dp_plots))
            ggsave(file=paste0(model$name, "/dp_", names(dp_plots)[i], ".pdf"), dp_plots[[i]])
    } else if(n_chains == 0){
        ## perfect simulation - start at all in one or all in different block configs
        ## all in one block
        init_sbm <- sbm(blocks(rep(1, Edges$numnodes)), rparams(1, model$params))
        lower    <- sampler(Edges, model, n_steps, "dp", sigma, currsbm=init_sbm)
        ## all in seperate blocks
        init_sbm <- sbm(blocks(1:Edges$numnodes), rparams(Edges$numnodes, model$params))
        upper    <- sampler(Edges, model, n_steps, "dp", sigma, currsbm=init_sbm)
        ## plots
        perfectDP <- data.frame(x=1:n_steps, upper=upper$postk, lower=lower$postk)
        save(perfectDP, file=paste0(model$name, "/dp_perfect.RDa"))
        ggsave(paste0(model$name, "/dp_perfect_trace.pdf"),
               ggplot(perfectDP, aes(x, upper, lower)) +
               geom_line(aes(y = lower, colour = "lower")) +
               geom_line(aes(y = upper, colour = "upper")) +
               xlab("Iteration") + ylab("Number of blocks") +
               labs(color="Chain") +
               scale_y_log10(),
               width=5, height=5)
    } else{
        DP <- mclapply(1:n_chains, function(n_chains) sampler(Edges, model, 2*n_steps, "dp", sigma), mc.cores=n_chains)
        x <- matrix(NA, n_chains, n_steps)
        y <- matrix(NA, n_chains, n_steps)
        thetaMeans <- list()
        thetaVars <- list()
        for(j in 1:n_chains){
            for(i in 1:n_steps){
                if(!is.null(DP[[j]])){
                    x[j,i] <- mean(DP[[j]]$postt[1,,n_steps+i], na.rm=TRUE)
                    y[j,i] <- var(DP[[j]]$postt[1,, n_steps+i], na.rm=TRUE)
                }
                thetaMeans[[j]] <- as.mcmc(x[j,])
                thetaVars[[j]] <- as.mcmc(y[j,])
            }
        }
        DPmeans <- gelman.diag(thetaMeans)
        DPvars  <- gelman.diag(thetaVars)
        save(DPmeans, DPvars, file=paste0(model$name, "/dp_rubin_gelman.RDa"))
    }
}
