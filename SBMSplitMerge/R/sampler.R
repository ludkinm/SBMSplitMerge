#' @title top level sampler function
#' @param Edges edges object
#' @param Mod model object
#' @param nSteps number of steps to run mcmc for
#' @param algorithm choice of algorithm options are: "conjugate", "gibbs", "dp", "rj"
#' @param sigma random walk parameter for theta
#' @param statusfreq frequency to print some status info out - in number of iterations
#' @param currsbm initial state for sbm object - if NULL one is drawn from mod
#' @param ... additional parameters to pass to sbm for currsbm
#' @return postz traces for block assignments z
#' @return postt traces for theta
#' @return postk traces for number of blocks kappa
#' @return postn traces for number of occupied blocks
#' @return nsteps how many steps did the algo use
#' @return algorithm which algo was used
#' @export
sampler <- function(Edges, Mod, nSteps=1000, algorithm="rj",
                    sigma=0.5, statusfreq, currsbm, ...){
    if(missing(statusfreq))
        statusfreq <- nSteps

    step <- switch( algorithm ,
                   "conjugate" = sampler.conj,
                   "gibbs" = sampler.gibbs,
                   "dp" = sampler.dp,
                   "rj" = sampler.rj
                   )

    if(algorithm=="conjugate"){
        if(!is.function(Mod$marglike))
            stop("Please set a marginal likelihood function in Mod")
        warning("Selected conjugate Model - ensure the edge-state and parameter Models are conjugate")
    }

    if( algorithm == "gibbs" & !Mod$blocks$fixkappa)
        stop("Can't use the gibbs sampler unless block Model has a fixed kappa")

    N <- Edges$numnodes
    if(missing(currsbm))
        currsbm <- rsbm(N, Mod)
    dtheta <- currsbm$params$dtheta
    postz <- array(NA, c(N, nSteps))
    postt <- array(NA, c(dtheta, N+1, nSteps))
    postn    <- array(NA, c(N, nSteps))
    posttbar <- rep(NA, nSteps)
    posttvar <- rep(NA, nSteps)
    postl <- rep(NA, nSteps)
    postk <- rep(NA, nSteps)

    ## store sample 1
    postz[,1]   <- currsbm$blocks$z
    postt[,1,1] <- currsbm$params$theta0
    postt[,1+1:currsbm$kappa,1] <- t(currsbm$params$thetak)
    postk[1] <- currsbm$kappa
    postl[1] <- loglike(Edges, currsbm, Mod)
    postn[1:currsbm$kappa,1] <- currsbm$blocks$sizes

    for(s in 2:nSteps){
        currsbm <- step(currsbm, Edges, Mod, sigma, ...)
        ## store
        postk[s] <- currsbm$kappa
        postl[s] <- loglike(Edges, currsbm, Mod)
        postz[,s]   <- currsbm$blocks$z
        postt[,1,s] <- currsbm$params$theta0
        postt[,1+1:currsbm$kappa,s] <- t(currsbm$params$thetak)
        postn[1:currsbm$kappa,s] <- currsbm$blocks$sizes
        status(s, statusfreq)
    }

    list(postz=postz, postl=postl, postt=postt, postk=postk, postn=postn, nsteps=nSteps, algorithm=algorithm)
}

#' print the status every \code{n} steps
#' @param s step number
#' @param n frequency
status <- function(s, n=10)
    if( s %% n == 0)
        cat(s, "\n")

#' accept propsbm with the acceptance probability alpha
#' @param currsbm current sbm state
#' @param propsbm proposed sbm state
#' @param Edges edges object
#' @param Mod model object
#' @param logjac log jacobian of transformation of variables
#' @param logu log density for auxiliary variables
#' @param ... additional arguments to pass to loglike
#' @return new currsbm
accept <- function(currsbm, propsbm, Edges, Mod, logjac=0, logu=0, ...){
    logl <- loglike(Edges, propsbm, Mod, ...) - loglike(Edges, currsbm, Mod, ...)
    logp <- dparams(propsbm, Mod) - dparams(currsbm, Mod)
    logb <- dblocks(propsbm, Mod) - dblocks(currsbm, Mod)
    A <- exp(logl+logp+logb+logjac+logu)
    ## potential that A up like exp(Inf - Inf?)
    if(!is.nan(A) & (stats::runif(1) < A))
        currsbm <- propsbm
    currsbm
}

#' conjugate model sampler
#' @param currsbm the current state of the sampler
#' @param Edges edges object
#' @param Mod model list
sampler.conj <- function(currsbm, Edges, Mod){
    for(i in 1:currsbm$numnodes){
        znoi <- zmat(currsbm$blocks)[,-i,drop=FALSE]
        ei   <- Edges$Edges[i,-i]
        p <- Mod$marglike(znoi, ei, Mod$params)
        q <- condprior(currsbm$blocks, Mod$blocks, i)
        p <- normaliselogs(p)
        newblock <- rcat(1,p)
        currsbm <- updateblock.sbm(currsbm, newblock, i)
    }
    currsbm
}

#' Gibbs sampling for node assignments
#' @param currsbm the current state of the sampler
#' @param Edges edges object
#' @param Mod model list
#' @param sigma random walk parameter for theta
sampler.gibbs <- function(currsbm, Edges, Mod, sigma){
    currsbm <- drawblocks.gibbs(currsbm, Edges, Mod)
    currsbm <- drawparams(currsbm, Edges, Mod, sigma)
    currsbm
}

#' dirichlet process sampler
#' @param currsbm the current state of the sampler
#' @param Edges edges object
#' @param Mod model list
#' @param sigma random walk parameter for theta
sampler.dp <- function(currsbm, Edges, Mod, sigma){
    currsbm <- drawblocks.dp(currsbm, Edges, Mod)
    currsbm <- drawparams(currsbm, Edges, Mod, sigma)
    currsbm
}

#' reversible jump MCMC split-merge sampler
#' @param currsbm the current state of the sampler
#' @param Edges edges object
#' @param Mod model list
#' @param sigma random walk parameter for theta
#' @param rho propensity to add a block
sampler.rj <- function(currsbm, Edges, Mod, sigma, rho=10){
    currsbm <- drawparams(currsbm, Edges, Mod, sigma=sigma)
    if( stats::runif(1) < 0.5 | currsbm$kappa == 1 ){
        currsbm <- splitavg(currsbm, Edges, Mod)
    } else{
        currsbm <- mergeavg(currsbm, Edges, Mod)
    }
    emptyblocks <- sum(currsbm$blocks$sizes==0)
    pdel <- emptyblocks/(emptyblocks+rho)
    if( stats::runif(1) < pdel ){
        currsbm <- delblock(currsbm, Edges, Mod, rho=rho)
    } else{
        currsbm <- addblock(currsbm, Edges, Mod, rho=rho)
    }
    currsbm <- drawblocks.gibbs(currsbm, Edges, Mod)
    currsbm
}
