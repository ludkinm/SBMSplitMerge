#' remove a block
#' @param z factor vector of block memberships
#' @param k a block (level of z)
#' @return blocks object with no block k
rmblock <- function(z, k)
    ## Takes a vector of block assignments and removes label k.
    ## Returned as block object
    blocks(factor(z, levels(z)[-k]))

#' add a block move
#' @param sbm an sbm object
#' @param edges an edges object
#' @param mod a model object
#' @param rho probability of choosing to add a block
#' @return an updated sbm object
addblock <- function(sbm, edges, mod, rho=1){
    ## proposes adding an empty block labelled kappa+1 to the sbm object
    ## auto accept-rejects
    ## returns sbm object
    nempty <- sum(sbm$blocks$sizes == 0)
    newparam <- rparams(1, mod$params)$thetak
    p <- params(sbm$params$theta0, rbind(sbm$params$thetak, newparam))
    b <- blocks(sbm$blocks$z, sbm$kappa+1)
    propsbm <- sbm(b,p)
    logb <- dblocks(propsbm, mod) - dblocks(sbm, mod)
    logu <- log(rho + nempty) - log(rho) - log(rho + nempty + 1)
    A <- exp(logb+logu)
    ## potential that A up like exp(Inf - Inf?)
    if(!is.nan(A) & (runif(1) < A))
        sbm <- propsbm
    sbm
}

#' delete a block move
#' @param sbm an sbm object
#' @param edges an edges object
#' @param mod a model object
#' @param rho probability of choosing to add a block
#' @return an updated sbm object
delblock <- function(Sbm, Edges, Mod, rho=1){
    ## proposes deleting an empty block (chosen at random among empty blocks)
    ## auto accept-rejects
    ## returns sbm object
    emptyind <- Sbm$blocks$sizes == 0
    nempty   <- sum(emptyind)
    k        <- which(emptyind)[sample(nempty,1)]
    b <- rmblock(Sbm$blocks$z, k)
    p <- params(Sbm$params$theta0, Sbm$params$thetak[-k, ,drop=FALSE])
    propsbm <- sbm(b,p)
    logb <- dblocks(propsbm, Mod) - dblocks(Sbm, Mod)
    logu <- log(rho + nempty) + log(rho) - log(rho + nempty - 1)
    A <- exp(logb+logu)
    ## potential that A up like exp(Inf - Inf?)
    if(!is.nan(A) & (runif(1) < A))
        Sbm <- propsbm
    Sbm
}

#' merge move using average to merge parameters
#' @param sbm an sbm object
#' @param edges an edges object
#' @param mod a model object
#' @param ... additional parameter to 'accept'
#' @return an updated sbm object
mergeavg <- function(sbm, edges, mod, ...){
    ## merge two blocks together chosen at random
    ## returns sbm
    kl <- rcat(2, rep(1, sbm$kappa), FALSE)
    k <- min(kl)
    l <- max(kl)
    mp <- mergeparams(sbm$params, k, l, mod$params)
    mb <- mergeblocks(sbm, k, l, edges, mod)
    propsbm <- sbm(mb$prop, mp$prop)
    logu <- log(sbm$kappa) - log(2) + (sbm$kappa==2)*log(2)
    accept(sbm, propsbm, edges, mod, mp$loga + mb$loga, logu, ...)
}

#' split move using average to merge parameters
#' @param sbm an sbm object
#' @param edges an edges object
#' @param mod a model object
#' @param ... additional parameter to 'accept'
#' @return an updated sbm object
splitavg <- function(Sbm, Edges, Mod, ...){
    ## split a block chosen at random
    ## returns sbm
    k       <- sample(Sbm$kappa, 1)
    ## propose parameters
    sp      <- splitparams(Sbm$params, k, Mod$params)
    pparams <- sp$prop
    sb      <- splitblocks(Sbm$blocks, pparams, k, Edges, Mod)
    pblocks <- sb$prop
    propsbm <- sbm(pblocks, pparams)
    ## acceptence probability
    logu <- log(2) - log(Sbm$kappa+1)  - (Sbm$kappa==1)*log(2)
    accept(Sbm, propsbm, Edges, Mod, logjac=sp$loga + sb$loga, logu=logu, ...)
}

#' merge move block merging
#' @param sbm an sbm object
#' @param k,l blocks to merge
#' @param edges an edges object
#' @param mod a model object
#' @param ... additional parameter to 'accept'
#' @return list(proposed block structure, log-acceptance-prob)
mergeblocks <- function(sbm, k, l, edges, mod, ...){
    ## merges two blocks, whilst calculating the inverse probability
    origz <- currz <- propz <- sbm$blocks$z

    ## make proposed block structure by relabelling any l's to k's
    ind   <- (origz == k) | (origz == l)
    propz[ind] <- k
    propblocks <- rmblock(propz, l)

    ## calculate reverse probability of assignment
    nodes <- which(ind)[sample(sum(ind))] # shuffle nodes
    currz[ind] <- NA
    currblocks <- blocks(currz, sbm$kappa) # running block model
    loga <- 0
    for(i in nodes){
        p <- nodelike(edges, currblocks, sbm$params, mod, i)
        p <- normaliselogs(p[c(k,l)])
        loga <- loga + log(p[origz[i] == c(k,l)])
        currz[i] <- origz[i]
    }
    list(prop = propblocks, loga = loga)
}

#' split move block mergingsplitting
#' @param Blocks a blocks object
#' @param pparams propoed_params
#' @param k block to split
#' @param Edges an edges object
#' @param mod a model object
#' @param ... additional parameter to 'accept'
#' @return list(proposed block structure, log-acceptance-prob)
splitblocks <- function(Blocks, pparams, k, Edges, Mod, ...){
    ## splits two block calculating the acceptance probability along the way
    ind  <- Blocks$z == k
    loga <- 0
    l <- Blocks$kappa+1
    if( sum(ind) > 0 ){                 # if k is not-empty
        propz <- Blocks$z
        propz[ind] <- NA
        ## calculate reverse probability of assignment
        nodes <- which(ind)[sample(sum(ind))] # shuffle nodes
        propblocks <- blocks(propz, l) # running block model
        for(i in nodes){
            p <- nodelike(Edges, propblocks, pparams, Mod, i)
            p <- normaliselogs(p[c(k,l)])
            j <- rcat(1,p)
            propblocks <- updateblock(propblocks, c(k,l)[j], i)
            loga <- loga - log(p[j])
        }
        ## make the larger block stay labelled as k
        if( propblocks$sizes[l] > propblocks$sizes[k] ){
            propz <- propblocks$z
            in_l <- propz == l
            in_k <- propz == k
            propz[in_l] <- k
            propz[in_k] <- l
            propblocks <- blocks(propz, l)
        }
    } else {
        propblocks <- blocks(Blocks$z, l)
    }
    list(prop = propblocks, loga = loga)
}


#' split move split params
splitparams <- function(x, ...)
    UseMethod("splitparams", x)

#' split move split params
#' @param params a parameter object to split
#' @param k block to split
#' @param pm parammod object
#' @return list(proposed_params, log-acceptance-prob)
splitparams.params <- function(params, k, pm){
    ## split params, draw x (w in paper) uniform and u normally distributed
    x    <- rbeta(params$dtheta, 5, 5)
    u    <- rnorm(params$dtheta, 0, 1)
    loga <- sum(dnorm(u, 0, 1, log=T))
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
#' @param pm parammod object
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
mergeparams <- function(x,...)
    UseMethod("mergeparams", x)

#' Merge step - parameter merging
#' @param params a params object
#' @param k,l blocks to merge
#' @param pm a parammod object
#' @return list(proposed_params, log-acceptance-prob)
mergeparams.default <- function(params, k, l, pm){
    ## gi
    x <- runif(params$dtheta)
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
#' @param pm a parammod object
#' @return list(proposed_params, log-acceptance-prob)
mergeparams.numeric <- function(thetak, thetal, x, pm){
    ## given thetak thetal and x (w in paper) and the parameter model
    ## return the merged value theta with log jacobian value
    theta <- pm$invt(x*pm$t(thetak) + (1-x)*pm$t(thetal))
    loga <- -pm$loggradt(theta) + pm$loggradt(thetak) + pm$loggradt(thetal) + sum(log(x)) + sum(log(1-x)) + length(x) * log(2)
    list(prop = theta, loga=loga)
}
