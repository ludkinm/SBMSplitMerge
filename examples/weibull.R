library(SBMSplitMerge)
edge_model <- edgemod(
        function(e, p) dweibull(e, p[1,,], p[2,,], log=TRUE)
        ,
        function(p) rweibull(1, p[1], p[2])
)
block_model <- dma(1, 10)
param_model <- parammod(
        function(params){
                dgamma(params$theta0[1], 1, 1, log=TRUE) +
                        dgamma(params$theta0[2], 1, 1, log=TRUE) +
                        sum(dgamma(params$thetak[,1], 1, 1, log=TRUE)) +
                        sum(dgamma(params$thetak[,2], 1, 1, log=TRUE))
        },
        function(kappa){
                params(
                        c(rgamma(1, 1, 1), rgamma(1, 1, 1))
                ,
                        cbind(rgamma(kappa, 1, 1), rgamma(kappa, 1, 1))
                )
        },
        function(x){ cbind(log(x[1]), log(x[2]))},
        function(x){ cbind(exp(x[1]), exp(x[2]))},
        function(x){ -log(x[1])-log(x[2])}
        )
model <- sbmmod(block_model, param_model, edge_model)
set.seed(1)
true_blocks   <- blocks(rep(c(1, 2, 3), c(10, 20, 20)))
true_params   <- params(c(1, 0.5), cbind(1, c(3,4,5)))
true_sbm      <- sbm(true_blocks, true_params)
weibull_edges <- redges(true_sbm, model$edge)
print(true_blocks)
plot(true_blocks)
plot(weibull_edges)

## run split merge
set.seed(1)
sm_output <- sampler(weibull_edges, model, 200, "rj", 0.5, 25)
sm_plots <- eval_plots(sm_output)
sm_plots$blocks_trace
sm_plots$post_pairs

## run DP
model <- sbmmod(crp(10), model$param, model$edge)
set.seed(1)
dp_output <- sampler(weibull_edges, model, 200, "dp", 0.5, 25)
dp_plots <- eval_plots(dp_output)
dp_plots$blocks_trace
dp_plots$post_pairs

## Run gibbs
model <- sbmmod(multinom(1, 3), model$param, model$edge)
set.seed(1)
gibbs_output <- sampler(edges=weibull_edges, sbmmod=model, nSteps=200, algorithm = "gibbs", sigma=0.1, 25)
gibbs_plots <- eval_plots(gibbs_output)
gibbs_plots$blocks_trace
gibbs_plots$post_pairs

MAPsbm <- function(output){
    map_index <- which.max(output$postl)
    kappa <- output$postk[map_index]
    sbm(
        blocks(output$postz[,map_index], kappa)
       ,
        params(output$postt[,1,map_index], t(output$postt[,1 + 1:kappa,map_index]))
    )
}

smMAP <- MAPsbm(sm_output)
dpMAP <- MAPsbm(dp_output)
gibbsMAP <- MAPsbm(gibbs_output)

plot(weibull_edges, smMAP)
plot(weibull_edges, dpMAP)
plot(weibull_edges, gibbsMAP)
