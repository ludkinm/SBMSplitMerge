#' remove a block
#' @param z factor vector of block memberships
#' @param k a block (level of z)
#' @return Blocks object with no block k
rmblock <- function(z, k)
    ## Takes a vector of block assignments and removes label k.
    ## Returned as block object
    blocks(factor(z, levels(z)[-k]))

#' add a block move
#' @param Sbm an SBM object
#' @param Edges an Edges object
#' @param Mod a model object
#' @param rho probability of choosing to add a block
#' @return an updated SBM object
addblock <- function(Sbm, Edges, Mod, rho=1){
    ## proposes adding an empty block labelled kappa+1 to the SBM object
    ## auto accept-rejects
    ## returns Sbm object
    nempty <- sum(Sbm$blocks$sizes == 0)
    newparam <- rparams(1, Mod$params)$thetak
    p <- params(Sbm$params$theta0, rbind(Sbm$params$thetak, newparam))
    b <- blocks(Sbm$blocks$z, Sbm$kappa+1)
    propSbm <- sbm(b,p)
    logb <- dblocks(propSbm, Mod) - dblocks(Sbm, Mod)
    logu <- log(rho + nempty) - log(rho) - log(rho + nempty + 1)
    A <- exp(logb+logu)
    ## potential that A up like exp(Inf - Inf?)
    if(!is.nan(A) & (stats::runif(1) < A))
        Sbm <- propSbm
    Sbm
}

#' delete a block move
#' @param Sbm an Sbm object
#' @param Edges an Edges object
#' @param Mod a model list
#' @param rho probability of choosing to add a block
#' @return an updated Sbm object
delblock <- function(Sbm, Edges, Mod, rho=1){
    ## proposes deleting an empty block (chosen at random among empty Blocks)
    ## auto accept-rejects
    ## returns Sbm object
    emptyind <- Sbm$blocks$sizes == 0
    nempty   <- sum(emptyind)
    k        <- which(emptyind)[sample(nempty,1)]
    b <- rmblock(Sbm$blocks$z, k)
    p <- params(Sbm$params$theta0, Sbm$params$thetak[-k, ,drop=FALSE])
    propSbm <- sbm(b,p)
    logb <- dblocks(propSbm, Mod) - dblocks(Sbm, Mod)
    logu <- log(rho + nempty) + log(rho) - log(rho + nempty - 1)
    A <- exp(logb+logu)
    ## potential that A up like exp(Inf - Inf?)
    if(!is.nan(A) & (stats::runif(1) < A))
        Sbm <- propSbm
    Sbm
}

#' merge move using average to merge parameters
#' @param Sbm an Sbm object
#' @param Edges an Edges object
#' @param Mod a model object
#' @param ... additional parameter to 'accept'
#' @return an updated Sbm object
mergeavg <- function(Sbm, Edges, Mod, ...){
    ## merge two Blocks together chosen at random
    ## returns Sbm
    kl <- rcat(2, rep(1, Sbm$kappa), FALSE)
    k <- min(kl)
    l <- max(kl)
    mp <- mergeparams(Sbm$params, k, l, Mod$params)
    mb <- mergeblocks(Sbm, k, l, Edges, Mod)
    propSbm <- sbm(mb$prop, mp$prop)
    logu <- log(Sbm$kappa) - log(2) + (Sbm$kappa==2)*log(2)
    accept(Sbm, propSbm, Edges, Mod, mp$loga + mb$loga, logu, ...)
}

#' split move using average to merge parameters
#' @param Sbm an Sbm object
#' @param Edges an Edges object
#' @param Mod a model object
#' @param ... additional parameter to 'accept'
#' @return an updated Sbm object
splitavg <- function(Sbm, Edges, Mod, ...){
    ## split a block chosen at random
    ## returns Sbm
    k       <- sample(Sbm$kappa, 1)
    ## propose parameters
    sp      <- splitparams(Sbm$params, k, Mod$params)
    pparams <- sp$prop
    sb      <- splitblocks(Sbm$blocks, pparams, k, Edges, Mod)
    pblocks <- sb$prop
    propSbm <- sbm(pblocks, pparams)
    ## acceptence probability
    logu <- log(2) - log(Sbm$kappa+1)  - (Sbm$kappa==1)*log(2)
    accept(Sbm, propSbm, Edges, Mod, logjac=sp$loga + sb$loga, logu=logu, ...)
}

#' merge move block merging
#' @param Sbm an Sbm object
#' @param k Blocks to merge
#' @param l Blocks to merge
#' @param Edges an Edges object
#' @param Mod a model object
#' @return list(proposed block structure, log-acceptance-prob)
mergeblocks <- function(Sbm, k, l, Edges, Mod){
    ## merges two Blocks, whilst calculating the inverse probability
    origz <- currz <- propz <- Sbm$blocks$z

    ## make proposed block structure by relabelling any l's to k's
    ind   <- (origz == k) | (origz == l)
    propz[ind] <- k
    propBlocks <- rmblock(propz, l)

    ## calculate reverse probability of assignment
    nodes <- which(ind)[sample(sum(ind))] # shuffle nodes
    currz[ind] <- NA
    currBlocks <- blocks(currz, Sbm$kappa) # running block model
    loga <- 0
    for(i in nodes){
        p <- nodelike(Edges, currBlocks, Sbm$params, Mod, i)
        p <- normaliselogs(p[c(k,l)])
        loga <- loga + log(p[origz[i] == c(k,l)])
        currz[i] <- origz[i]
    }
    list(prop = propBlocks, loga = loga)
}

#' split move block mergingsplitting
#' @param Blocks a Blocks object
#' @param pparams propoed_params
#' @param k block to split
#' @param Edges an Edges object
#' @param Mod a model object/list
#' @return list(proposed block structure, log-acceptance-prob)
splitblocks <- function(Blocks, pparams, k, Edges, Mod){
    ## splits two block calculating the acceptance probability along the way
    ind  <- Blocks$z == k
    loga <- 0
    l <- Blocks$kappa+1
    if( sum(ind) > 0 ){                 # if k is not-empty
        propz <- Blocks$z
        propz[ind] <- NA
        ## calculate reverse probability of assignment
        nodes <- which(ind)[sample(sum(ind))] # shuffle nodes
        propBlocks <- blocks(propz, l) # running block model
        for(i in nodes){
            p <- nodelike(Edges, propBlocks, pparams, Mod, i)
            p <- normaliselogs(p[c(k,l)])
            j <- rcat(1,p)
            propBlocks <- updateblock(propBlocks, c(k,l)[j], i)
            loga <- loga - log(p[j])
        }
        ## make the larger block stay labelled as k
        if( propBlocks$sizes[l] > propBlocks$sizes[k] ){
            propz <- propBlocks$z
            in_l <- propz == l
            in_k <- propz == k
            propz[in_l] <- k
            propz[in_k] <- l
            propBlocks <- blocks(propz, l)
        }
    } else {
        propBlocks <- blocks(Blocks$z, l)
    }
    list(prop = propBlocks, loga = loga)
}


#' split move split params
#' @param x object for dispatch
#' @param ... additional arguments for method
splitparams <- function(x, ...)
    UseMethod("splitparams", x)

#' split move split params
#' @param params a parameter object to split
#' @param k block to split
#' @param pm paramMod object
#' @return list(proposed_params, log-acceptance-prob)
splitparams.params <- function(params, k, pm){
    ## split params, draw x (w in paper) uniform and u normally distributed
    x    <- stats::rbeta(params$dtheta, 5, 5)
    u    <- stats::rnorm(params$dtheta, 0, 1)
    loga <- sum(stats::dnorm(u, 0, 1, log=T))
    l    <- params$kappa+1
    propthetak <- rbind(params$thetak,0)
    theta <- params$thetak[k,]
    tmp   <- splitparams(theta, u, x, pm)
    loga  <- loga + tmp$loga
    propthetak[c(k,l),] <- tmp$prop
    propp <- params(params$theta0, propthetak)
    list(prop = propp, loga=loga)
}

#' split move split params
#' @param theta a parameter to split
#' @param u auxiliary variable
#' @param x auxiliary variable
#' @param pm paramMod object
#' @return list(proposed_params, log-acceptance-prob)
splitparams.numeric <- function(theta, u, x, pm){
    ## given value for theta u and x and parameter model - perform the split to yield theta_k and theta_l
    ## returns list with proposed param structure and log-jacobian of transformation
    propk <- pm$invt((pm$t(theta) + u)/2/x)
    propl <- pm$invt((pm$t(theta) - u)/2/(1-x))
    loga <- pm$loggradt(theta) - pm$loggradt(propk) - pm$loggradt(propl) - sum(log(x)) - sum(log(1-x)) - length(x) * log(2)
    prop <- rbind(propk, propl)
    list(prop = prop, loga=loga)
}

#' merge parameters
#' @param x an object to dispatch on
#' @param ... additional arguments for methods
mergeparams <- function(x,...)
    UseMethod("mergeparams", x)

#' Merge step - parameter merging
#' @param params a params object
#' @param k Blocks to merge
#' @param l Blocks to merge
#' @param pm a paramMod object
#' @return list(proposed_params, log-acceptance-prob)
mergeparams.default <- function(params, k, l, pm){
    ## gi
    x <- stats::runif(params$dtheta)
    propthetak <- params$thetak
    tmp <- mergeparams(params$thetak[k,], params$thetak[l,], x, pm)
    propthetak[k,] <- tmp$prop
    loga <- tmp$loga
    propthetak <- propthetak[-l,,drop=FALSE]
    propp <- params(params$theta0, propthetak)
    list(prop = propp, loga = loga)
}

#' Merge step - parameter merging
#' @param thetak,thetal parameters to merge
#' @param x auxiliary parameter
#' @param pm a paramMod object
#' @return list(proposed_params, log-acceptance-prob)
mergeparams.numeric <- function(thetak, thetal, x, pm){
    ## given thetak thetal and x (w in paper) and the parameter model
    ## return the merged value theta with log jacobian value
    theta <- pm$invt(x*pm$t(thetak) + (1-x)*pm$t(thetal))
    loga <- -pm$loggradt(theta) + pm$loggradt(thetak) + pm$loggradt(thetal) + sum(log(x)) + sum(log(1-x)) + length(x) * log(2)
    list(prop = theta, loga=loga)
}
