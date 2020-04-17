#' Class for edge data
#' @param e a matrix or array representing the raw edge-state data
#' @param sym is the network symetric? (is edge-state ji = edgestate ji)
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
#' @param ... additional arguments to append to edgemod internal list
#' @return an edgemod object
#' @export
edgemod <- function(dg, rg, ...){
    out <- list(dg=dg, ...)
    if(!missing(rg))
        out$rg <- rg
    class(out) <- append(class(out), "edgemod")
    out
}

#' print an edges object
#' @param x an edges object
#' @param ... (unused)
#' @export
print.edges <- function(x, ...)
    print(paste("edges object on", x$numnodes, "nodes"))

#' check if an object is an edges object
#' @param x an R object
#' @return TRUE if x is an edges object
#' @export
is.edges <- function(x)
    inherits(x,"edges")

#' plots edges objects
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot .data geom_raster theme scale_alpha xlab ylab aes element_blank
#' @param x an edges object
#' @param Blocks a blocks object or sbm object
#' @param sorted sort by block membership in blocks before plotting?
#' @param xlab label for x axis
#' @param ylab label for y axis
#' @param ... parameters for image
#' @return ggplot2 plot of edges in a raster
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
#' @param x an R object
#' @param ... additional arguments
#' @export
redges <- function(x, ...)
    UseMethod("redges", x)

#' simulate edges
#' @param x an sbm object
#' @param mod a model object
#' @param ... additional arguments passed to redges.params
#' @return an edges object
#' @export
redges.sbm <- function(x, mod, ...)
    redges.params(x$params, x$blocks, mod$edges,...)

#' simulate edges
#' @param x a params object
#' @param blocks a blocks object
#' @param edgemod an edgemod object
#' @param sym should the network be symmetric?
#' @param loops should the network have self-loops?
#' @param ... additional arguments (unused)
#' @return an edges object
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
#' @param edges an edges object
#' @param x an R object for dispatch
#' @param ... additional arguments
#' @return matrix same size as edges$E with density of each edge
#' @export
dedges <- function(edges, x, ...)
    UseMethod("dedges", x)

#' density of edges
#' @param edges an edges object
#' @param x an sbm object
#' @param mod a model object
#' @param ... additional arguments for dedges.params
#' @return matrix same size as edges$E with density of each edge
dedges.sbm <- function(edges, x, mod, ...)
    dedges.params(edges, x$params, x$blocks, mod$edges, ...)

#' density of edges
#' @param edges an edges object
#' @param x an params object
#' @param blocks a blocks object
#' @param edgemod an edgemod object
#' @param ... additional arguments passed to dedges.numeric
#' @return matrix same size as edges$E with density of each edge
#' @export
dedges.params <- function(edges, x, blocks, edgemod, ...){
    pmat <- parammat(blocks, x)
    dedges(edges, pmat, edgemod, ...)
}

#' likelihood of edges
#' @param edges an edges object
#' @param x a matrix of parameters
#' @param edgemod an edgemod object
#' @param ... additional arguments passed to edgemod$dg
#' @return likelihood of edges under the edgemod using parameters in matrix pmat which has same size as edgse$E
#' @export
dedges.numeric <- function(edges, x, edgemod, ...)
    edgemod$dg(edges$E, x, ...)

#' likelihood of edges
#' @param edges an edges object
#' @param x an object to dispatch on
#' @param ... additional arguments for inherited method
#' @return numeric likelihood of edges under x
#' @export
loglike <- function(edges, x, ...)
    UseMethod("loglike", x)

#' likelihood of edges
#' @param edges an edges object
#' @param x an sbm object
#' @param mod a model object
#' @param ... additional arguments for loglike.params
#' @return numeric likelihood of edges under x
#' @export
loglike.sbm <- function(edges, x, mod, ...)
    loglike.params(edges, x$params, x$blocks, mod$edges, ...)

#' likelihood of edges
#' @param edges an edges object
#' @param x a params object
#' @param blocks a blocks object
#' @param edgemod an edgemod object
#' @param na.rm passed to sum to remove NAs or not (default FALSE)
#' @param ... additional arguments for dedges
#' @return likelihood of edges under the edgemod using parameters in matrix pmat which has same size as edgse$E
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

#' get the edges eminating from node i
#' @param x object for dispatch
#' @param ... additional arguments
edgesfromi <- function(x,...)
    UseMethod("edgesfromi", x)

edgesfromi.dynedges <- function(edges, i, kappa=1)
    dynedges(edges$E[,rep(i,kappa), -i, drop=FALSE], edges$obtimes)

#' get the edges eminating from node i
#' @param Edges an edges object
#' @param i the node incident to the returned edges
#' @param kappa number of blocks
#' @return an edges object containing only the edges incident to i
edgesfromi.edges <- function(Edges, i, kappa=1)
    edges(Edges$E[rep(i,kappa), -i, drop=FALSE], Edges$sym, Edges$loops)


#' get likelihood of edges eminating from node i
#' @param Edges an edges object
#' @param x an object
#' @param ... additional arguments
#' @return likelihood of edges eminating from node i
#' @export
nodelike <- function(Edges, x, ...)
    UseMethod("nodelike", x)

#' get likelihood of edges eminating from node i
#' @param Edges an edges object
#' @param x an sbm object
#' @param Mod a model object
#' @param i the node of interest
#' @param ... additional arguments for nodelike.blocks
#' @return likelihood of edges eminating from node i
#' @export
nodelike.sbm <- function(Edges, x, Mod, i, ...)
    nodelike.blocks(Edges, x$blocks, x$params, Mod, i, ...)

#' get likelihood of edges eminating from node i
#' @param Edges an edges object
#' @param x an blocks object
#' @param Params a params object
#' @param Mod a model object
#' @param i the node of interest
#' @param ... additional arguments passed to dedges
#' @return likelihood of edges eminating from node i
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


#' Make an edgemod model with bernoulli edge-states
#' @param ... additional parameters to pass to rbinom
#' @return an edgemod representing edgestates which are Bernoulli distributed
#' @export
edges_bern <- function(...)
    edgemod(
        ...
       ,
        rg = function(p, ...) stats::rbinom(1, 1, p, ...)
       ,
        dg = function(e, p) stats::dbinom(e, 1, p, log=TRUE)
    )

#' Make an edgemod model with Poisson edge-states
#' @param ... additional parameters to pass to rpois
#' @return an edgemod representing edgestates which are Poisson distributed
#' @export
edges_pois <- function(...)
    edgemod(
        ...
       ,
        rg = function(p, ...) stats::rpois(1, p, ...)
       ,
        dg = function(e, p) stats::dpois(e, p, log=TRUE)
    )

#' Make an edgemod model with Normal edge-states
#' @param ... additional parameters to pass to rnorm
#' @return an edgemod representing edgestates which are Gaussian distributed
#' @export
edges_norm <- function(...)
    edgemod(
        ...
       ,
        rg = function(p, ...) stats::rnorm(1, p[1], p[2], ...)
       ,
        dg = function(e, p) stats::dnorm(e, p[1,,], p[2,,], log=TRUE)
    )

#' Make an edgemod model with NegativeBinomial edge-states
#' @param ... additional parameters to pass to rnbinom
#' @return an edgemod representing edgestates which are NegativeBinomial distributed
#' @export
edges_nbin <- function(...)
    edgemod(
        ...
       ,
        rg = function(p, ...) stats::rnbinom(1, p[1], p[2], ...)
       ,
        dg = function(e, p) stats::dnbinom(e, p[1,,], p[2,,], log=TRUE)
    )
