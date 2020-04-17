#' simulate parameters for the given model in a Metropolis scheme
#' @param currsbm current \code{sbm} object
#' @param Edges edge state data \code{edges} object
#' @param Mod model object for the \code{sbm}
#' @param sigma parameter for \code{drawparam}
#' @return updated \code{sbm} object
#' @export
drawparams <- function(currsbm, Edges, Mod, sigma=0.1){
    for(k in 0:currsbm$kappa){
        tmp <- drawparam(currsbm, k, Mod, sigma)
        currsbm <- accept(currsbm, tmp$prop, Edges, Mod, tmp$logjac)
    }
    currsbm
}

#' draw parameters
#' @param x object for dispatch
#' @param ... additional arguments for method
#' @export
drawparam <- function(x, ...)
    UseMethod("drawparam", x)

#' simulate parameter for the kth parameter
#' @param x current \code{sbm} object
#' @param k index of parameter to update
#' @param Mod model object for the \code{sbm} in \code{x}
#' @param sigma parameter for \code{drawparam}
#' @param ... additional arguments for \code{drawparam.params}
#' @return \code{list(proposed sbm, log-jacobian term)}
#' @export
drawparam.sbm <- function(x, k, Mod, sigma, ...){
    p <- drawparam(x$params, x$blocks, k, Mod$params, sigma, ...)
    list(prop=sbm(x$blocks, p$prop), logjac=p$logjac)
}

#' simulate parameter for the kth parameter
#' @param x current \code{params} object
#' @param Blocks current \code{blocks} object
#' @param k index of parameter to update
#' @param Parammod model object for the parameters
#' @param sigma parameter for \code{rw}
#' @param ... additional arguments (unused)
#' @return \code{list(proposed params, log-jacobian term)}
#' @export
drawparam.params <- function(x, Blocks, k, Parammod, sigma, ...){
    if(k == 0){
        ## if operating on the between-block parameter...
        if(Blocks$kappa > 1){
            ## ...and the number of blocks is more than 1
            ## do a random walk on the between-block parameter
            tmp <- rw(x$theta0, Parammod, sigma)
        } else{
            ## ...and the number of blocks is 1
            ## draw from the prior
            tmp <- rparam(0, Parammod)
        }
        x$theta0 <- tmp$prop
    } else{
        ## if operating on a within block parameter...
        if(Blocks$sizes[k] > 1){
            ## ...and the block has more than 1 node
            ## do a random walk on the parameter
            tmp <- rw(x$thetak[k,], Parammod, sigma)
        } else{
            ## else draw from the prior
            tmp <- rparam(k, Parammod)
        }
        x$thetak[k,] <- tmp$prop
    }
    list(prop=x, logjac = tmp$logjac)
}

#' random walk
#' @param p a parameter
#' @param pm a \code{parammod} object
#' @param sigma - scale of random walk
#' @return \code{ist(proposed parameter, locjacobian)}
rw <- function(p, pm, sigma){
    pprime <- pm$invt(pm$t(p) + stats::rnorm(length(p), 0, sigma))
    logjac <- pm$loggradt(pprime) - pm$loggradt(p)
    list(prop = pprime, logjac=logjac)
}

#' draw from the prior
#' @param k index of parameter
#' @param pm a \code{parammod} object
#' #' @return \code{ist(proposed parameter, locjacobian=0)}
rparam <- function(k, pm){
    if(k==0){
        prop <- rparams(0, pm)$theta0
    } else{
        prop <- rparams(1, pm)$thetak[1,]
    }
    list(prop=prop, logjac=0)
}
