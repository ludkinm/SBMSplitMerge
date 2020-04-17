#' Class for edge data
#' @param e a matrix or array representing the raw edge-state data
#' @param sym is the network symmetric? (\code{e[ji] = e[ji]})
#' @param loops does the network contain self-loops? (edges from node i to i)
#' @param ... additional arguments to append to edges internal list
#' @return an edges object
#' @export
edges <- function(e, sym, loops, ...){
    if(missing(sym))
        sym <- all(e==t(e))
    if(missing(loops))
        loops <- !all(diag(e)==0)
    out <- list(
        E = e
       ,
        numnodes = ncol(e)
       ,
        sym = sym
       ,
        loops = loops
      ,
        ...
    )
    class(out) <- append("edges", class(out))
    out
}

#' Class for edge models
#' @param dg function to calculate likelihood of edges
#' @param rg function to simulate edges optional
#' @param ... additional arguments to append to \code{edgemod} internal list
#' @return an \code{edgemod} object
#' @export
edgemod <- function(dg, rg, ...){
    out <- list(dg=dg, ...)
    if(!missing(rg))
        out$rg <- rg
    class(out) <- append(class(out), "edgemod")
    out
}

#' print an \code{edges} object
#' @param x an \code{edges} object
#' @param ... (unused)
#' @export
print.edges <- function(x, ...)
    print(paste("edges object on", x$numnodes, "nodes"))

#' check if an object is an \code{edges} object
#' @param x an R object
#' @return TRUE if x is an \code{edges} object
#' @export
is.edges <- function(x)
    inherits(x,"edges")

#' plots edges objects
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot .data geom_raster theme scale_alpha xlab ylab aes element_blank
#' @param x an \code{edges} object
#' @param Blocks a blocks object or \code{sbm} object
#' @param sorted sort by block membership in \code{blocks} before plotting?
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param ... parameters for \code{image}
#' @return \code{ggplot2} plot of edges in a raster
#' @export
plot.edges <- function(x, Blocks, sorted=TRUE, xlab="Node", ylab="Node", ...){
    ord <- 1:x$numnodes
    if(sorted & !missing(Blocks)){
        if(is.sbm(Blocks))
            Blocks <- Blocks$blocks
        ord <- order(Blocks$z)
        Blocks$z <- Blocks$z[ord]
    }
    colnames(x$E) <- rownames(x$E) <- NULL
    df <- reshape2::melt(x$E[ord,ord])
    if(!missing(Blocks)){
        df$block1 <- 0
        df$block2 <- 0
        for(i in 1:x$numnodes)
            df$block2[df$Var2 == i] <- df$block1[df$Var1 == i] <- as.numeric(as.character(Blocks$z[i]))
        df$block <- df$block1 * (df$block1 == df$block2)
    }
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$Var1, .data$Var2)) +
        ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
              panel.grid.major=ggplot2::element_blank()) +
        ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
    if(!missing(Blocks)){
        p <- p + ggplot2::geom_raster(ggplot2::aes(fill=factor(.data$block)))
    }
    p <- p + ggplot2::geom_raster(ggplot2::aes(alpha=as.numeric(.data$value))) +
        ggplot2::scale_alpha(range=c(0.4,0), name="value")
    p
}

#' @rdname plot.edges
#' @export
image.edges <- plot.edges

#' convert x to a symmetric matrix by taking the lower triangle
#' @param x object for dispatch
#' @param ... additional arguments for method
makesymmetric <- function(x, ...)
    UseMethod("makesymmetric", x)

makesymmetric.matrix <- function(e){
    e[lower.tri(e)] <- t(e)[lower.tri(e)]
    e
}

makesymmetric.array <- function(e){
    for(i in 1:dim(e)[1])
        e[i,,] <- makesymmetric(e[i,,])
    e
}

#' generate an edges object from x
#' @param x an R object for dispatch
#' @param ... additional arguments
#' @export
redges <- function(x, ...)
    UseMethod("redges", x)

#' simulate edges
#' @param x an \code{sbm} object
#' @param mod a model list
#' @param ... additional arguments passed to \code{redges.params}
#' @return an edges object
#' @export
redges.sbm <- function(x, mod, ...)
    redges.params(x$params, x$blocks, mod$edges,...)

#' simulate edges
#' @param x a \code{params} object
#' @param blocks a \code{blocks} object
#' @param edgemod an \code{edgemod} object
#' @param sym should the network be symmetric?
#' @param loops should the network have self-loops?
#' @param ... additional arguments (unused)
#' @return an \code{edges} object
#' @export
redges.params <- function(x, blocks, edgemod, sym=T, loops=F, ...){
    if(is.null(edgemod$rg))
        stop("The edgemodel does not specify a random method via $rg so can't simulate edges")
    pmat <- parammat(blocks, x)
    x <- apply(pmat, 2:3, edgemod$rg)
    if(sym)
        x <- makesymmetric(x)
    if(!loops)
        diag(x) <- 0
    edges(x)
}

#' density of edges
#' @param edges an \code{edges} object
#' @param x an R object for dispatch
#' @param ... additional arguments
#' @return matrix same size as \code{edges$E} with density of each edge
#' @export
dedges <- function(edges, x, ...)
    UseMethod("dedges", x)

#' density of edges
#' @param edges an \code{edges} object
#' @param x an \code{sbm} object
#' @param mod a model list
#' @param ... additional arguments for \code{dedges.params}
#' @return matrix same size as \code{edges$E} with density of each edge
dedges.sbm <- function(edges, x, mod, ...)
    dedges.params(edges, x$params, x$blocks, mod$edges, ...)

#' density of edges
#' @param edges an \code{edges} object
#' @param x an \code{params} object
#' @param blocks a \code{blocks} object
#' @param edgemod an \code{edgemod} object
#' @param ... additional arguments passed to \code{dedges.numeric}
#' @return matrix same size as \code{edges$E} with density of each edge
#' @export
dedges.params <- function(edges, x, blocks, edgemod, ...){
    pmat <- parammat(blocks, x)
    dedges(edges, pmat, edgemod, ...)
}

#' likelihood of edges
#' @param edges an \code{edges} object
#' @param x a matrix of parameters (with same size as \code{edgse$E})
#' @param edgemod an \code{edgemod} object
#' @param ... additional arguments passed to \code{edgemod$dg}
#' @return likelihood of edges under the \code{edgemod} using parameters in matrix \code{pmat}
#' @export
dedges.numeric <- function(edges, x, edgemod, ...)
    edgemod$dg(edges$E, x, ...)

#' likelihood of edges
#' @param edges an \code{edges} object
#' @param x an object to dispatch on
#' @param ... additional arguments for inherited method
#' @return numeric likelihood of edges under \code{x}
#' @export
loglike <- function(edges, x, ...)
    UseMethod("loglike", x)

#' likelihood of edges
#' @param edges an \code{edges} object
#' @param x an \code{sbm} object
#' @param mod a model list
#' @param ... additional arguments for \code{loglike.params}
#' @return numeric likelihood of edges under \code{x}
#' @export
loglike.sbm <- function(edges, x, mod, ...)
    loglike.params(edges, x$params, x$blocks, mod$edges, ...)

#' likelihood of edges
#' @param edges an \code{edges} object
#' @param x a \code{params} object
#' @param blocks a \code{blocks} object
#' @param edgemod an \code{edgemod} object
#' @param na.rm passed to sum to remove \code{NA}s or not (default \code{FALSE})
#' @param ... additional arguments for \code{dedges}
#' @return likelihood of \code{edges} under the \code{edgemod} using parameters in \code{x}
#' @export
loglike.params <- function(edges, x, blocks, edgemod, na.rm=FALSE, ...){
    x <- dedges(edges, x, blocks, edgemod, ...)
    if(!edges$loops)
        diag(x) <- 0
    ll <- sum(x, na.rm=na.rm)
    if(edges$sym)
        ll <- ll/2
    ll
}

#' get the edges emanating from node i
#' @param x object for dispatch
#' @param ... additional arguments
edgesfromi <- function(x,...)
    UseMethod("edgesfromi", x)

edgesfromi.dynedges <- function(edges, i, kappa=1)
    dynedges(edges$E[,rep(i,kappa), -i, drop=FALSE], edges$obtimes)

#' get the edges emanating from node i
#' @param Edges an \code{edges} object
#' @param i the node incident to the returned \code{edges}
#' @param kappa number of blocks
#' @return an \code{edges} object containing only the edges incident to i
edgesfromi.edges <- function(Edges, i, kappa=1)
    edges(Edges$E[rep(i,kappa), -i, drop=FALSE], Edges$sym, Edges$loops)


#' get likelihood of edges emanating from node i
#' @param Edges an \code{edges} object
#' @param x an object
#' @param ... additional arguments
#' @return likelihood of edges emanating from node i
#' @export
nodelike <- function(Edges, x, ...)
    UseMethod("nodelike", x)

#' get likelihood of edges emanating from node i
#' @param Edges an \code{edges} object
#' @param x an \code{sbm} object
#' @param Mod a \code{model} object
#' @param i the node of interest
#' @param ... additional arguments for \code{nodelike.blocks}
#' @return likelihood of edges emanating from node i
#' @export
nodelike.sbm <- function(Edges, x, Mod, i, ...)
    nodelike.blocks(Edges, x$blocks, x$params, Mod, i, ...)

#' get likelihood of edges emanating from node i
#' @param Edges an \code{edges} object
#' @param x an \code{blocks} object
#' @param Params a \code{params} object
#' @param Mod a model list
#' @param i the node of interest
#' @param ... additional arguments passed to \code{dedges}
#' @return likelihood of edges emanating from node i
#' @export
nodelike.blocks <- function(Edges, x, Params, Mod, i, ...){
    Edgemod <- Mod$edges
    kap <- x$kappa + !Mod$blocks$fixkappa
    znoti <- zmat(x, kap)[,-i, drop=FALSE]
    ei <- edgesfromi(Edges, i, kap)
    zi <- diag(kap)
    p  <- rep(0, kap)
    pmat <- parammat(zi, znoti, Params)
    rowSums(dedges(ei, pmat, Edgemod, ...), na.rm=TRUE)
}


#' @title Bernoulli edge model
#' @description Make an \code{edgemod} model with Bernoulli edge-states
#' @param ... additional parameters to pass to \code{rbinom}
#' @return an \code{edgemod}
#' @export
edges_bern <- function(...)
    edgemod(
        ...
       ,
        rg = function(p, ...) stats::rbinom(1, 1, p, ...)
       ,
        dg = function(e, p) stats::dbinom(e, 1, p, log=TRUE)
    )

#' @title Poisson edge model
#' @description Make an \code{edgemod} model with Poisson edge-states
#' @param ... additional parameters to pass to \code{rpois}
#' @return an \code{edgemod}
#' @export
edges_pois <- function(...)
    edgemod(
        ...
       ,
        rg = function(p, ...) stats::rpois(1, p, ...)
       ,
        dg = function(e, p) stats::dpois(e, p, log=TRUE)
    )

#' @title Normal edge model
#' @description Make an \code{edgemod} model with Normal edge-states
#' @param ... additional parameters to pass to \code{rnorm}
#' @return an \code{edgemod}
#' @export
edges_norm <- function(...)
    edgemod(
        ...
       ,
        rg = function(p, ...) stats::rnorm(1, p[1], p[2], ...)
       ,
        dg = function(e, p) stats::dnorm(e, p[1,,], p[2,,], log=TRUE)
    )

#' @title Negative-Binomial edge model
#' @description Make an \code{edgemod} model with Negative-Binomial edge-states
#' @param ... additional parameters to pass to \code{rnbinom}
#' @return an \code{edgemod}
#' @export
edges_nbin <- function(...)
    edgemod(
        ...
       ,
        rg = function(p, ...) stats::rnbinom(1, p[1], p[2], ...)
       ,
        dg = function(e, p) stats::dnbinom(e, p[1,,], p[2,,], log=TRUE)
    )
