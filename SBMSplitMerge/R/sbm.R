#' Class \code{sbm}
#' @param blocks a \code{blocks} object
#' @param params a \code{params} object
#' @return an \code{sbm} object
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

#' print an \code{sbm} object
#' @param x an \code{sbm} object
#' @param ... additional arguments for formatting
#' @export
print.sbm <- function(x,...)
    cat(format(x,...), "\n")

format.sbm <- function(sbm,...)
    c("SBM object:\n\n",format(sbm$blocks,...), "\n\n", format(sbm$params,...))

#' plot an \code{sbm} object as an \code{image}
#' @param x an \code{sbm} object
#' @param col colours as in \code{plot.default}
#' @param axes axes as in \code{plot.default}
#' @param ... additional arguments for plot
#' @seealso plot.default
#' @export
plot.sbm <- image.sbm <- function(x, col, axes=FALSE, ...){
    if(missing(col))
        col <- c(0,rainbow(x$kappa))
    plot(x$blocks, axes=axes, col=col, ...)
    plot(x$params, col=rep(col, each=x$params$dtheta), ...)
}

#' check if an R object is an \code{sbm} object
#' @param x an R object
#' @return TRUE if x is an \code{sbm} object
#' @export
is.sbm <- function(x)
    inherits(x, "sbm")

#' simulate an \code{sbm}
#' @param numnodes number of nodes
#' @param mod a model list
#' @return an \code{sbm} object
#' @export
rsbm <- function(numnodes, mod)
    sbm(tmp <- rblocks(numnodes, mod$blocks), rparams(tmp$kappa, mod$params))
