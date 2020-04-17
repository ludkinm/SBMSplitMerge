#' calculate the marginal likelihood for a node for samplers using conjugate models
#' @title Marginal likelihood model for Bernoulli distributed edges
#' @param znoi a matrix of blockassignments without node i
#' @param ei edge-states incident to i
#' @param parammod a parammod object
#' @export
marglike_bern <- function(znoi, ei, parammod){
    ## block sizes
    nk  <- rowSums(znoi)
    mk <- c(ei %*% t(znoi))
    mk <- rbind(mk, nk - mk)
    if(!any(nk == 0))
        mk <- cbind(mk, 0)
    m0 <- c(sum(ei),sum(1-ei)) - mk
    a0 <- parammod$alpha + m0
    ak <- parammod$beta + mk
    logp <- lbeta(ak[1,],ak[2,]) + lbeta(a0[1,],a0[2,])
    logp
}

#' calculate the marginal likelihood for a node for samplers using conjugate models
#' @title Marginal likelihood model for Poisson distributed edges
#' @param znoi a matrix of blockassignments without node i
#' @param ei edge-states incident to i
#' @param parammod a parammod object
#' @export
marglike_pois <- function(znoi, ei, parammod){
    ## block sizes
    nk <- rowSums(znoi)
    mk <- rbind(ei %*% t(znoi), nk)
    if(!any(nk == 0))
        mk <- cbind(mk, 0)
    m0 <- rowSums(mk) - mk
    a0 <- parammod$alpha + m0
    ak <- parammod$beta + mk
    logp <- lgamma(ak[1,]) - ak[1,] * log(ak[2,]) + lgamma(a0[1,]) - a0[1,] * log(a0[2,])
    logp
}

#' calculate the marginal likelihood for a node for samplers using conjugate models
#' @title Marginal likelihood model for Normal distributed edges
#' @param znoi a matrix of blockassignments without node i
#' @param ei edge-states incident to i
#' @param parammod a parammod object
#' @export
marglike_norm <- function(znoi, ei, parammod){
    rho0 <- parammod$alpha[1]
    nu0  <- parammod$alpha[2]
    be0  <- parammod$alpha[3]
    al0  <- parammod$alpha[4]
    rhok <- parammod$beta[1]
    nuk  <- parammod$beta[2]
    bek  <- parammod$beta[3]
    alk  <- parammod$beta[4]
    nk   <- rowSums(znoi)
    mk   <- c(ei %*% t(znoi))
    ssk  <- c(ei^2 %*% t(znoi))
    if(!any(nk == 0)){
        nk <- c(nk, 0)
        mk <- c(mk, 0)
        ssk <- c(ssk, 0)
    }
    m0 <- sum(mk) - mk
    n0 <- sum(nk) - nk
    ss0 <- sum(ssk) - ssk
    postbek <- bek + 0.5*ssk - 0.5*mk^2/nk + nuk*(mk/nk - rhok)^2/(nuk+nk)/nk/2
    postbe0 <- be0 + 0.5*ss0 - 0.5*m0^2/n0 + nu0*(m0/n0 - rho0)^2/(nu0+n0)/n0/2
    postbek[is.nan(postbek)] <- bek
    postbe0[is.nan(postbe0)] <- be0
    logp <- 0
    logp <- logp + lgamma(alk + nk/2) - lgamma(alk)
    logp <- logp + 0.5*log(nuk) - 0.5*log(nuk + nk)
    logp <- logp - nk/2 * log(2*pi)
    logp <- logp + alk*log(bek)
    logp <- logp - (alk + nk/2)*log(postbek)
    logp <- logp + lgamma(al0 + n0/2) - lgamma(al0)
    logp <- logp + 0.5*log(nu0) - 0.5*log(nu0 + n0)
    logp <- logp - n0/2 * log(2*pi)
    logp <- logp + al0*log(be0)
    logp <- logp - (al0 + n0/2)*log(postbe0)
    logp
}
