#' Class sbm
#' @param blocks a blocks object
#' @param params a params object
#' @return an sbm object
#' @export
sbm <- function(blocks, params){
    if(blocks$kappa != params$kappa)
        stop("Mis-matched number of blocks")
    rownames(params$thetak) <- names(blocks$n)
    out <- list(
        blocks = blocks
       ,
        params = params
       ,
        kappa = blocks$kappa
       ,
        numnodes = blocks$numnodes
    )
    class(out) <- c(class(out), "sbm")
    out
}

#' print an SBM object
#' @param x an SBM object
#' @param ... additional arguments for formatting
#' @export
print.sbm <- function(x,...)
    cat(format(x,...), "\n")

format.sbm <- function(sbm,...)
    c("SBM object:\n\n",format(sbm$blocks,...), "\n\n", format(sbm$params,...))

#' plot an SBM object as an image
#' @param x an SBM object
#' @param col colors as in plot.default
#' @param axes axes as in plot.default
#' @param ... additional arguments to plot
#' @seealso plot.default
#' @export
plot.sbm <- image.sbm <- function(x, col, axes=FALSE, ...){
    if(missing(col))
        col <- c(0,rainbow(x$kappa))
    plot(x$blocks, axes=axes, col=col, ...)
    plot(x$params, col=rep(col, each=x$params$dtheta), ...)
}

#' check if an R object is an SBM object
#' @param x an R object
#' @return TRUE if x is an SBM object
#' @export
is.sbm <- function(x)
    inherits(x, "sbm")

#' simulate an sbm
#' @param numnodes umber of nodes
#' @param mod a model object
#' @return an sbm object
#' @export
rsbm <- function(numnodes, mod)
    sbm(tmp <- rblocks(numnodes, mod$blocks), rparams(tmp$kappa, mod$params))
