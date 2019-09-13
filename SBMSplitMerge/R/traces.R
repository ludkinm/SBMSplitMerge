
#' helper function for trace plots
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_raster theme scale_x_continuous scale_y_continuous element_blank
postimage <- function(mat, discrete=FALSE){
    df <- reshape2::melt(mat)
    a <- ggplot2::aes(Var2, Var1, fill=value)
    if(discrete)
        a <- ggplot2::aes(Var2, Var1, fill=factor(value))
    p <- ggplot2::ggplot(df, a) + ggplot2::geom_raster() +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
                       panel.grid.major=ggplot2::element_blank())
    p
}

#' plot a trace of the blocks from MCMC samples
#' @importFrom ggplot2 xlab ylab guides guide_legend
#' @param postz output from sampler
#' @return ggplot object of node assignments against iteration number
#' @export
blocktrace <- function(postz){
    p <- postimage(postz, TRUE)
    p <- p + ggplot2::ylab("Node") + ggplot2::xlab("Step") + ggplot2::guides(fill=ggplot2::guide_legend(title="Block"))
    p
}

#' plot a trace of the number of blocks from MCMC samples
#' @importFrom ggplot2 ggplot aes geom_line
#' @param postk output from sampler
#' @return ggplot object of kappa against iteration number
#' @export
numblockstrace <- function(postk){
    ggplot2::ggplot(data=data.frame(step=seq_along(postk), kappa=postk), ggplot2::aes(x=step, y=kappa)) + ggplot2::geom_line()
}


#' plot a trace of parameter values from MCMC samples
#' @importFrom reshape2 melt
#' @importFrom ggplot2 xlab ylab guides guide_legend geom_hline
#' @param theta output from sampler
#' @param truetheta if supplied, adds guidelines for the true parameters
#' @return ggplot object of theta value against iteration number
#' @export
paramtrace <- function(theta, truetheta){
    ind <- colSums(apply(theta, 2, is.na)) != dim(theta)[3]
    theta[,ind,][is.na(theta[,ind,])] <- 0
    dft <- reshape2::melt(theta[,ind,,drop=FALSE], id=Dimension)
    names(dft) <- c("Dimension", "Block", "Step", "Value")
    p <- ggplot2::ggplot(dft, ggplot2::aes(x=Step, y=Value, colour=factor(Block))) + ggplot2::geom_line()
    if(!missing(truetheta))
        p <- p + ggplot2::geom_hline(yintercept=truetheta)
    p
}

#' mean proportion of times two nodes were in the same block under MCMC samples
#' @param postz output from sampler
#' @return matrix P with P[i,j] = proportion of times i and j are in the same block under postz
postpairs <- function(postz){
    N <- nrow(postz)
    P <- matrix(0,N,N)
    for(i in 2:N)
        for(j in 1:(i-1))
            P[j,i] <- P[i,j] <- mean(postz[i,] == postz[j,])
    P
}

#' modal block assignments from MCMC samples
#' @export
#' @param postz output from sampler
#' @return a blocks object with the modal block assignments under postz
modeblocks <- function(postz)
    blocks(apply(postz, 1, function(x) which.max(tabulate(x, max(postz)))))


#' get a set of evaluation plots from MCMC samples
#' @export
#' @param output from sampler
#' @param burnin burn-in period (a vector of iteration numbers to subset outputs)
#' @return list of ggplot objects (with descriptive names)
eval_plots <- function(output, burnin, truetheta){
    if(missing(burnin))
        burnin <- 1:output$nsteps
    nb <- numblockstrace(output$postk[burnin])
    pp <- postpairs(output$postz[,burnin])
    pp_plot <- postimage(pp)
    ## order for posterior plotting
    ind <- order(colSums(pp))
    pp_sorted <- postimage(pp[ind,ind])
    pt <- paramtrace(output$postt[,,burnin, drop=FALSE], truetheta)
    pt_sorted <- paramtrace(output$postt[,ind,burnin, drop=FALSE], truetheta)
    bt <- blocktrace(output$postz[,burnin])
    bt_sorted <- blocktrace(output$postz[ind,burnin])
    list(num_blocks_trace = nb, post_pairs = pp_plot, post_pairs_sorted = pp_sorted, blocks_trace = bt, param_trace = pt, param_trace_sorted = pt_sorted, blocks_trace_sorted = bt_sorted, sortind = ind, pp=pp)
}
