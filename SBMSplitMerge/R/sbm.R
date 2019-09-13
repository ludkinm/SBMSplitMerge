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

#' @export
print.sbm <- function(sbm,...)
    cat(format(sbm,...), "\n")

#' @export
format.sbm <- function(sbm)
    c("SBM object:\n\n",format(sbm$blocks), "\n\n", format(sbm$params))

#' @export
plot.sbm <- image.sbm <- function(sbm, col, axes=FALSE,...){
    if(missing(col))
        col <- c(0,rainbow(sbm$kappa))
    plot(sbm$blocks, axes=axes, col=col, ...)
    plot(sbm$params, col=rep(col, each=sbm$params$dtheta), ...)
}

#' @export
is.sbm <- function(sbm)
    inherits(sbm, "sbm")

#' simulate an sbm
#' @param numnode number of nodes
#' @param mod a model object
#' @return an sbm object
#' @export
rsbm <- function(numnodes, mod, ...)
    sbm(tmp <- rblocks(numnodes, mod$blocks, ...), rparams(tmp$kappa, mod$params))
