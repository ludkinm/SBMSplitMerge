#' Class for edge data
#' @param e a matrix or array representing the raw edge-state data
#' @param sym is the network symetric? (is edge-state ji = edgestate ji)
#' @param loops does the network contain self-loops? (edges from node i to i)
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
#' @param rg function to simulate edges
#' @param dg function to calculate likelihood of edges
#' @return an edgemod object
#' @export
edgemod <- function(rg, dg, ...){
    out <- list(
        rg=rg
       ,
        dg=dg
    )
    class(out) <- append(class(out), "edgemod")
    out
}

#' @export
print.edges <- function(edges)
    print(paste("edges object on", edges$numnodes, "nodes"))

#' @export
is.edges <- function(edges)
    inherits(edges,"edges")

#' plots edges objects
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_raster theme scale_fill_manual scale_alpha xlab ylab aes element_blank
#' @param edges an edges object
#' @param blocks a blocks object or sbm object
#' @param sorted sort by block membership in blocks before plotting?
#' @param ... parameters for image
#' @export
plot.edges <- function(edges, blocks, sorted=TRUE, xlab="Node", ylab=" ", ...){
    ord <- 1:edges$numnodes
    if(sorted & !missing(blocks)){
        if(is.sbm(blocks))
            blocks <- blocks$blocks
        ord <- order(blocks$z)
        blocks$z <- blocks$z[ord]
    }
    colnames(edges$E) <- rownames(edges$E) <- NULL
    df <- reshape2::melt(edges$E[ord,ord])
    if(!missing(blocks)){
        df$block1 <- 0
        df$block2 <- 0
        for(i in 1:edges$numnodes)
            df$block2[df$Var2 == i] <- df$block1[df$Var1 == i] <- as.numeric(as.character(blocks$z[i]))
        df$block <- df$block1*(df$block1 == df$block2)
    }
    p <- ggplot2::ggplot(df, ggplot2::aes(Var1, Var2)) +
        ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
              panel.grid.major=ggplot2::element_blank()) +
        ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
    if(!missing(blocks)){
        p <- p + ggplot2::geom_raster(ggplot2::aes(fill=factor(block)))
    }
    p <- p + ggplot2::geom_raster(ggplot2::aes(alpha=as.numeric(value))) +
        ggplot2::scale_alpha(range=c(0.4,0), name="value")
    p
}

#' @export
image.edges <- plot.edges

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
#' @export
redges <- function(x, ...)
    UseMethod("redges", x)

#' simulate edges
#' @param sbm an sbm object
#' @param mod a model object
#' @return an edges object
#' @export
redges.sbm <- function(sbm, mod, ...)
    redges.params(sbm$params, sbm$blocks, mod$edges,...)

#' simulate edges
#' @param params a params object
#' @param blocks a blocks object
#' @param edgemod an edgemod object
#' @param sym should the network be symmetric?
#' @param loops should the network have self-loops?
#' @return an edges object
#' @export
redges.params <- function(params, blocks, edgemod, sym=T, loops=F, ...){
    pmat <- parammat(blocks, params)
    x <- apply(pmat, 2:3, edgemod$rg)
    if(sym)
        x <- makesymmetric(x)
    if(!loops)
        diag(x) <- 0
    edges(x)
}

#' density of edges
#' @param edges an edges object
#' @param x an object
#' @return matrix same size as edges$E with density of each edge
#' @export
dedges <- function(edges, x, ...)
    UseMethod("dedges", x)

#' density of edges
#' @param edges an edges object
#' @param sbm an sbm object
#' @param mod a model object
#' @return matrix same size as edges$E with density of each edge
dedges.sbm <- function(edges, sbm, mod, ...)
    dedges.params(edges, sbm$params, sbm$blocks, mod$edges, ...)

#' density of edges
#' @param edges an edges object
#' @param params an params object
#' @param blocks a blocks object
#' @param edgemod an edgemod object
#' @return matrix same size as edges$E with density of each edge
#' @export
dedges.params <- function(edges, params, blocks, edgemod){
    pmat <- parammat(blocks, params)
    dedges(edges, pmat, edgemod)
}

#' likelihood of edges
#' @param edges an edges object
#' @param pmat a matrix of parameters
#' @param edgemod an edgemod object
#' @return likelihood of edges under the edgemod using parameters in matrix pmat which has same size as edgse$E
#' @export
dedges.numeric <- function(edges, pmat, edgemod)
    edgemod$dg(edges$E, pmat)

#' likelihood of edges
#' @param edges an edges object
#' @param x an object
#' @return numeric likelihood of edges under x
#' @export
loglike <- function(edges, x, ...)
    UseMethod("loglike", x)

#' likelihood of edges
#' @param edges an edges object
#' @param sbm an sbm object
#' @param mod a model object
#' @return numeric likelihood of edges under x
#' @export
loglike.sbm <- function(edges, sbm, mod, ...)
    loglike.params(edges, sbm$params, sbm$blocks, mod$edges)

#' likelihood of edges
#' @param edges an edges object
#' @param params a params object
#' @param blocks a blocks object
#' @param edgemod an edgemod object
#' @param na.rm=F passed to sum
#' @return likelihood of edges under the edgemod using parameters in matrix pmat which has same size as edgse$E
#' @export
loglike.params <- function(edges, params, blocks, edgemod, na.rm=F){
    x <- dedges(edges, params, blocks, edgemod)
    if(!edges$loops)
        diag(x) <- 0
    ll <- sum(x, na.rm=na.rm)
    if(edges$sym)
        ll <- ll/2
    ll
}

#' get the edges eminating from node i
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
#' @return likelihood of edges eminating from node i
#' @export
nodelike <- function(edges, x, ...)
    UseMethod("nodelike", x)

#' get likelihood of edges eminating from node i
#' @param Edges an edges object
#' @param sbm an sbm object
#' @param Mod a model object
#' @param i the node of interest
#' @return likelihood of edges eminating from node i
#' @export
nodelike.sbm <- function(Edges, Sbm, Mod, i, ...)
    nodelike.blocks(Edges, Sbm$blocks, Sbm$params, Mod, i, ...)

#' get likelihood of edges eminating from node i
#' @param Edges an edges object
#' @param blocks an blocks object
#' @param Params a params object
#' @param Mod a model object
#' @param i the node of interest
#' @return likelihood of edges eminating from node i
#' @export
nodelike.blocks <- function(Edges, Blocks, Params, Mod, i, ...){
    Edgemod <- Mod$edges
    kap <- Blocks$kappa + !Mod$blocks$fixkappa
    znoti <- zmat(Blocks, kap)[,-i, drop=FALSE]
    ei <- edgesfromi(Edges, i, kap)
    zi <- diag(kap)
    p  <- rep(0, kap)
    pmat <- parammat(zi, znoti, Params)
    rowSums(dedges(ei, pmat, Edgemod), na.rm=TRUE)
}


#' Make an edgemod model with bernoulli edge-states
#' @param ... additional parameters to pass to rbinom
#' @return an edgemod representing edgestates which are Bernoulli distributed
#' @export
edges_bern <- function(...)
    edgemod(
        ...
       ,
        rg = function(p, ...) rbinom(1, 1, p, ...)
       ,
        dg = function(e, p) dbinom(e, 1, p, log=TRUE)
    )

#' Make an edgemod model with Poisson edge-states
#' @param ... additional parameters to pass to rpois
#' @return an edgemod representing edgestates which are Poisson distributed
#' @export
edges_pois <- function(...)
    edgemod(
        ...
       ,
        rg = function(p, ...) rpois(1, p, ...)
       ,
        dg = function(e, p) dpois(e, p, log=TRUE)
    )

#' Make an edgemod model with Normal edge-states
#' @param ... additional parameters to pass to rnorm
#' @return an edgemod representing edgestates which are Gaussian distributed
#' @export
edges_norm <- function(...)
    edgemod(
        ...
       ,
        rg = function(p, ...) rnorm(1, p[1], p[2], ...)
       ,
        dg = function(e, p) dnorm(e, p[1,,], p[2,,], log=TRUE)
    )

#' Make an edgemod model with NegativeBinomial edge-states
#' @param ... additional parameters to pass to rnbinom
#' @return an edgemod representing edgestates which are NegativeBinomial distributed
#' @export
edges_nbin <- function(...)
    edgemod(
        ...
       ,
        rg = function(p, ...) rnbinom(1, p[1], p[2], ...)
       ,
        dg = function(e, p) dnbinom(e, p[1,,], p[2,,], log=TRUE)
    )
