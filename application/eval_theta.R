## table for latex
library(coda)
eval_theta <- function(output_theta, range=1:5, burnin){
    dimtheta <- dim(output_theta)[1]
    thetas <- output_theta[,range,burnin,drop=FALSE]
    qs <- apply(thetas, 1:2, quantile, probs=c(0.05, 0.5, 0.95))
    ess <- t(apply(thetas, 1:2, effectiveSize))
    tab <- NULL
    for(k in 1:dimtheta)
        tab <- cbind(tab, t(qs[,k,]), ESS=ess[,k])
    ps <- list()
    length(ps) <- dimtheta
    for(k in 1:dimtheta){
        df <- data.frame(burnin, cbind(t(thetas[k,,])))
        names(df) <- c("Index", paste0("Theta", 1:ncol(thetas)-1))
        mf <- reshape2::melt(df, id="Index")
        names(mf) <- c("Iteration", "var", "Value")
        ps[[k]] <- ggplot(mf, aes(x=Iteration, y=Value, col=var)) + geom_line()
    }
    list(theta_table=tab, plots=ps)
}
