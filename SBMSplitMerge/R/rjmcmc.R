#' remove a block
#' @param z factor vector of block memberships
#' @param k a block (level of z)
#' @return a \code{blocks} object without block k
rmblock <- function(z, k)
    ## Takes a vector of block assignments and removes label k.
    ## Returned as block object
    blocks(factor(z, levels(z)[-k]))

#' add a block move
#' @param SBM an \code{sbm} object
#' @param Edges an \code{edges} object
#' @param Mod a model list
#' @param rho probability of choosing to add a block
#' @return an updated \code{sbm} object
addblock <- function(SBM, Edges, Mod, rho=1){
    ## proposes adding an empty block labelled kappa+1 to the SBM object
    ## auto accept-rejects
    ## returns SBM object
    nempty <- sum(SBM$blocks$sizes == 0)
    newparam <- rparams(1, Mod$params)$thetak
    p <- params(SBM$params$theta0, rbind(SBM$params$thetak, newparam))
    b <- blocks(SBM$blocks$z, SBM$kappa+1)
    propSBM <- sbm(b,p)
    logb <- dblocks(propSBM, Mod) - dblocks(SBM, Mod)
    logu <- log(rho + nempty) - log(rho) - log(rho + nempty + 1)
    A <- exp(logb+logu)
    ## potential that A up like exp(Inf - Inf?)
    if(!is.nan(A) & (stats::runif(1) < A))
        SBM <- propSBM
    SBM
}

#' delete a block move
#' @param SBM an \code{sbm} object
#' @param Edges an \code{edges} object
#' @param Mod a model list
#' @param rho probability of choosing to add a block
#' @return an updated \code{sbm} object
delblock <- function(SBM, Edges, Mod, rho=1){
    ## proposes deleting an empty block (chosen at random among empty Blocks)
    ## auto accept-rejects
    ## returns SBM object
    emptyind <- SBM$blocks$sizes == 0
    nempty   <- sum(emptyind)
    k        <- which(emptyind)[sample(nempty,1)]
    b <- rmblock(SBM$blocks$z, k)
    p <- params(SBM$params$theta0, SBM$params$thetak[-k, ,drop=FALSE])
    propSBM <- sbm(b,p)
    logb <- dblocks(propSBM, Mod) - dblocks(SBM, Mod)
    logu <- log(rho + nempty) + log(rho) - log(rho + nempty - 1)
    A <- exp(logb+logu)
    ## potential that A up like exp(Inf - Inf?)
    if(!is.nan(A) & (stats::runif(1) < A))
        SBM <- propSBM
    SBM
}

#' merge move using average to merge parameters
#' @param SBM an \code{sbm} object
#' @param Edges an \code{edges} object
#' @param Mod a model list
#' @param ... additional parameter to 'accept'
#' @return an updated \code{sbm} object
mergeavg <- function(SBM, Edges, Mod, ...){
    ## merge two Blocks together chosen at random
    ## returns SBM
    kl <- rcat(2, rep(1, SBM$kappa), FALSE)
    k <- min(kl)
    l <- max(kl)
    mp <- mergeparams(SBM$params, k, l, Mod$params)
    mb <- mergeblocks(SBM, k, l, Edges, Mod)
    propSBM <- sbm(mb$prop, mp$prop)
    logu <- log(SBM$kappa) - log(2) + (SBM$kappa==2)*log(2)
    accept(SBM, propSBM, Edges, Mod, mp$loga + mb$loga, logu, ...)
}

#' split move using average to merge parameters
#' @param SBM an \code{sbm} object
#' @param Edges an \code{edges} object
#' @param Mod a model list
#' @param ... additional parameter to 'accept'
#' @return an updated \code{sbm} object
splitavg <- function(SBM, Edges, Mod, ...){
    ## split a block chosen at random
    ## returns SBM
    k       <- sample(SBM$kappa, 1)
    ## propose parameters
    sp      <- splitparams(SBM$params, k, Mod$params)
    pparams <- sp$prop
    sb      <- splitblocks(SBM$blocks, pparams, k, Edges, Mod)
    pblocks <- sb$prop
    propSBM <- sbm(pblocks, pparams)
    ## acceptence probability
    logu <- log(2) - log(SBM$kappa+1)  - (SBM$kappa==1)*log(2)
    accept(SBM, propSBM, Edges, Mod, logjac=sp$loga + sb$loga, logu=logu, ...)
}

#' merge move block merging
#' @param SBM an \code{sbm} object
#' @param k Blocks to merge
#' @param l Blocks to merge
#' @param Edges an \code{edges} object
#' @param Mod a model list
#' @return \code{list(proposed block structure, log-acceptance-prob)}
mergeblocks <- function(SBM, k, l, Edges, Mod){
    ## merges two Blocks, whilst calculating the inverse probability
    origz <- currz <- propz <- SBM$blocks$z

    ## make proposed block structure by relabelling any l's to k's
    ind   <- (origz == k) | (origz == l)
    propz[ind] <- k
    propBlocks <- rmblock(propz, l)

    ## calculate reverse probability of assignment
    nodes <- which(ind)[sample(sum(ind))] # shuffle nodes
    currz[ind] <- NA
    currBlocks <- blocks(currz, SBM$kappa) # running block model
    loga <- 0
    for(i in nodes){
        p <- nodelike(Edges, currBlocks, SBM$params, Mod, i)
        p <- normaliselogs(p[c(k,l)])
        loga <- loga + log(p[origz[i] == c(k,l)])
        currz[i] <- origz[i]
    }
    list(prop = propBlocks, loga = loga)
}

#' split move: blocks
#' @param Blocks a \code{blocks} object
#' @param pparams proposed \code{params }
#' @param k block to split
#' @param Edges an \code{edges} object
#' @param Mod a model list
#' @return \code{list(proposed block structure, log-acceptance-prob)}
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


#' split move: parameters
#' @param x object for dispatch
#' @param ... additional arguments for method
splitparams <- function(x, ...)
    UseMethod("splitparams", x)

#' split move: \code{params}
#' @param params a \code{params} object to split
#' @param k block to split
#' @param parammod \code{parammod} object
#' @return \code{list(proposed_params, log-acceptance-prob)}
splitparams.params <- function(params, k, parammod){
    ## split params, draw x (w in paper) uniform and u normally distributed
    x    <- stats::rbeta(params$dtheta, 5, 5)
    u    <- stats::rnorm(params$dtheta, 0, 1)
    loga <- sum(stats::dnorm(u, 0, 1, log=T))
    l    <- params$kappa+1
    propthetak <- rbind(params$thetak,0)
    theta <- params$thetak[k,]
    tmp   <- splitparams(theta, u, x, parammod)
    loga  <- loga + tmp$loga
    propthetak[c(k,l),] <- tmp$prop
    propp <- params(params$theta0, propthetak)
    list(prop = propp, loga=loga)
}

#' split move: \code{params}
#' @param theta a parameter to split
#' @param u auxiliary variable
#' @param x auxiliary variable
#' @param parammod \code{parammod} object
#' @return \code{list(proposed_params, log-acceptance-prob)}
splitparams.numeric <- function(theta, u, x, parammod){
    ## given value for theta u and x and parameter model - perform the split to yield theta_k and theta_l
    ## returns list with proposed param structure and log-jacobian of transformation
    propk <- parammod$invt((parammod$t(theta) + u)/2/x)
    propl <- parammod$invt((parammod$t(theta) - u)/2/(1-x))
    loga <- parammod$loggradt(theta) - parammod$loggradt(propk) - parammod$loggradt(propl) - sum(log(x)) - sum(log(1-x)) - length(x) * log(2)
    prop <- rbind(propk, propl)
    list(prop = prop, loga=loga)
}

#' merge parameters
#' @param x an object to dispatch on
#' @param ... additional arguments for methods
mergeparams <- function(x,...)
    UseMethod("mergeparams", x)

#' Merge step: parameters
#' @param params a \code{params object}
#' @param k Blocks to merge
#' @param l Blocks to merge
#' @param parammod a \code{parammod} object
#' @return \code{list(proposed_params, log-acceptance-prob)}
mergeparams.default <- function(params, k, l, parammod){
    ## gi
    x <- stats::runif(params$dtheta)
    propthetak <- params$thetak
    tmp <- mergeparams(params$thetak[k,], params$thetak[l,], x, parammod)
    propthetak[k,] <- tmp$prop
    loga <- tmp$loga
    propthetak <- propthetak[-l,,drop=FALSE]
    propp <- params(params$theta0, propthetak)
    list(prop = propp, loga = loga)
}

#' Merge step - parameter merging
#' @param thetak,thetal parameters to merge
#' @param x auxiliary parameter
#' @param parammod a \code{parammod} object
#' @return \code{list(proposed_params, log-acceptance-prob)}
mergeparams.numeric <- function(thetak, thetal, x, parammod){
    ## given thetak thetal and x (w in paper) and the parameter model
    ## return the merged value theta with log jacobian value
    theta <- parammod$invt(x*parammod$t(thetak) + (1-x)*parammod$t(thetal))
    loga <- -parammod$loggradt(theta) + parammod$loggradt(thetak) + parammod$loggradt(thetal) + sum(log(x)) + sum(log(1-x)) + length(x) * log(2)
    list(prop = theta, loga=loga)
}
