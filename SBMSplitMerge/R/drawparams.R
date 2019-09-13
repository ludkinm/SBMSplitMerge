#' simulate parameters for the given model in a Metropolis scheme
#' @param currsbm current state of sbm under the MCMC chain
#' @param Edges edge state data
#' @param Mod model object for the sbm
#' @param sigma parameter for drawparam
#' @return updated sbm object
#' @export
drawparams <- function(currsbm, Edges, Mod, sigma=0.1){
    for(k in 0:currsbm$kappa){
        tmp <- drawparam(currsbm, k, Mod, sigma)
        currsbm <- accept(currsbm, tmp$prop, Edges, Mod, tmp$logjac)
    }
    currsbm
}

#' @export
drawparam <- function(x, ...)
    UseMethod("drawparam", x)

#' simulate parameter for the kth parameter
#' @param currsbm current state of sbm under the MCMC chain
#' @param k index of parameter to update
#' @param Mod model object for the sbm
#' @param sigma parameter for drawparam
#' @return list(proposed sbm, log-jacobian term)
#' @export
drawparam.sbm <- function(currsbm, k, Mod, sigma){
    p <- drawparam(currsbm$params, currsbm$blocks, k, Mod$params, sigma)
    list(prop=sbm(currsbm$blocks, p$prop), logjac=p$logjac)
}

#' simulate parameter for the kth parameter
#' @param Params current params object
#' @param Blocks current blocks object
#' @param k index of parameter to update
#' @param Parammod model object for the parameters
#' @param sigma parameter for rw
#' @return list(proposed params, log-jacobian term)
#' @export
drawparam.params <- function(Params, Blocks, k, Parammod, sigma){
    propparams <- Params
    if(k == 0){
        if(Blocks$kappa > 1){
            tmp <- rw(Params$theta0, Parammod, sigma)
        } else{
            tmp <- rparam(0, Parammod)
        }
        propparams$theta0 <- tmp$prop
    } else{
        if(Blocks$sizes[k] > 1){
            tmp <- rw(Params$thetak[k,], Parammod, sigma)
        } else{
            tmp <- rparam(k, Parammod)
        }
        propparams$thetak[k,] <- tmp$prop
    }
    list(prop=propparams, logjac = tmp$logjac)
}

#' random walk
#' @param p a parameter
#' @param pm parammodel - contains the transform/map to project parameter to the real line
#' @param sigma - scale of random walk
rw <- function(p, pm, sigma){
    pprime <- pm$invt(pm$t(p) + rnorm(length(p), 0, sigma))
    logjac <- pm$loggradt(pprime) - pm$loggradt(p)
    list(prop = pprime, logjac=logjac)
}

#' draw from the prior - locjac = 0
#' @param k index of parameter
#' @param pm parammodel - contains the transform/map to project parameter to the real line and the prior on parameters
rparam <- function(k, pm){
    if(k==0){
        prop <- rparams(0, pm)$theta0
    } else{
        prop <- rparams(1, pm)$thetak[1,]
    }
    list(prop=prop, logjac=0)
}
