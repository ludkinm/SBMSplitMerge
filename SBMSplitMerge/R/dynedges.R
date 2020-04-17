## class
dynedges <- function(E, obtimes, sym){
    out <- list(
        E = E
       ,
        N = dim(E)[2]
       ,
        obtimes = obtimes
       ,
        sym = sym
       ,
        T = length(obtimes)
       ,
        dt = diff(obtimes)
    )
    class(out) <- append(class(out), "dynedges")
    out
}

## adjmatrix at a certain time
edgesattime <- function(dynedges, t)
    edges(dynedges$E[t,,], dynedges$sym)

## testers
is.dynedges <- function(dynedges)
    inherits(dynedges,"dynedges")

## printers
print.dynedges <- function(x,...)
    cat(format(x,...), "\n")

format.dynedges <- function(x, ...)
    c("Dynamic Edges object:\nN =", x$N, "\nT =", x$T, "\nobtimes:\n", x$obtimes)

## plotters
plot.dynedges <- image.dynedges <- function(dynedges, sbm, animate=F, ...){
    if(missing(sbm))
        sbm <- NULL
    plot(edgesattime(dynedges, 1), sbm, ...)
    if(animate)
        for(t in 1:length(dynedges$obtimes)){
            plot(edgesattime(dynedges, t), sbm, main=dynedges$obtimes[t],...)
            Sys.sleep(0.1)
        }
}

# simulate dynamic edges
rdynedges <- function(x, ...)
    UseMethod("rdynedges", x)

rdynedges.sbm <- function(sbm, obtimes, ...)
    rdynedges.params(sbm$params, sbm$blocks, obtimes, ...)

rdynedges.params <- function(params, blocks, obtimes, sym=TRUE, ...){
    x <- redges(params, parammat(blocks, params), blocks$N, obtimes, ...)
    if(sym)
        x <- makesymmetric(x)
    dynedges(x, obtimes, sym)
}

## ctsbm helpers
rctsbm <- function(obtimes, theta){
    ## single edge random function
    ## transition function
    phi <- theta[1]
    rho <- theta[2]
    trans <- function(now, state)
        now + stats::rexp(1, rates[state+1])
    nt    <- length(obtimes)
    rates  <- c(phi, 1-phi)*rho
    states <- stats::rbinom(1,1,phi) ## initial state
    jumptime <- obtimes[1]
    while(jumptime[1] < obtimes[nt]){
        ## then do a jump - ie swap states
        jumptime <- c(trans(jumptime[1], states[1]), jumptime)
        states   <- c(1-states[1], states)
    }
    stats::stepfun(rev(jumptime), rev(c(states,NA)))(obtimes)
}

dctsbm <- function(x, ...)
    UseMethod("dctsbm", x)

dctsbm.default <- function(eij, dt, theta)
    stats::dbinom(eij[1], 1, theta[1], log=T) + sum(stats::dbinom(eij[-1], 1, theta[1] + (eij[-length(eij)] - theta[1]) * exp(-dt * theta[2]), log=T))

dctsbm.matrix <- function(E, dt, theta)
    stats::dbinom(E[1,], 1, theta[1], log=T) + colSums(stats::dbinom(E[-1,], 1, theta[1] + (E[-dim(E)[1],] - theta[1]) * exp(-dt * theta[2]), log=T))

dctsbm.array <- function(E, dt, theta){
    last <- dim(E)[1]
    ind1 <- rep(1, last-1)
    ind2 <- rep(2, last-1)
    p <- theta[ind1,,] + (E[-last,,] - theta[ind1,,]) * exp(-dt * theta[ind2,,])
    stats::dbinom(E[1,,], 1, theta[1,,], log=T) + colSums(stats::dbinom(E[-1,,], 1, p, log=T))
}

dctsbm2 <- function(dt,Enow,Enext,phi,rho)
    stats::dbinom(Enext, 1, phi + (Enow - phi) * exp(-dt * rho), log=TRUE)
