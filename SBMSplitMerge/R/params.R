#' @title \code{parammod} S3 object
#' @description make a \code{parammod} object
#' @param alpha parameters for \code{theta_0}
#' @param beta parameters for \code{theta_k}
#' @param rf random function to draw parameters (the prior) \code{rf(n, x)} draws n parameters with hyperparameters x
#' @param df density function for parameters \code{df(theta, x)} density of theta with hyperparameters x
#' @param t map taking parameter space to real line
#' @param invt map taking real line to parameter space
#' @param loggradt log of the gradient of the map \code{t}
#' @return a \code{parammod} object
#' @export
parammod <- function(alpha, beta, rf, df, t, invt, loggradt){
    out <- list(alpha = alpha,
                beta = beta,
                rf = rf,
                df=df,
                t = t,
                invt = invt,
                loggradt = loggradt
                )
    class(out) <- unique(c(class(out), "parammod"))
    out
}

#' @title \code{params} S3 object
#' @description make a \code{params} object from the between-block parameter \code{theta0} and a vector of within block parameters \code{thetak}
#' @param theta0 between block parameter
#' @param thetak within block parameters
#' @return a \code{params} object
#' @export
params <- function(theta0, thetak){
    thetak <- as.matrix(thetak)
    theta0 <- as.matrix(theta0)
    out <- list(theta0 = theta0, thetak = thetak,
                kappa = nrow(thetak), dtheta = ncol(thetak))
    class(out) <- c(class(out), "params")
    out
}

#' @title generate parameters
#' @description Generate parameters from a give \code{parammod} object with \code{kappa} blocks
#' @param kappa number of blocks
#' @param parammod a \code{parammod} object
#' @return a \code{params} object
#' @export
rparams <- function(kappa, parammod){
    theta0 <- parammod$rf(1, parammod$alpha)
    thetak <- parammod$rf(kappa, parammod$beta)
    params(theta0, thetak)
}

#' @title Density of \code{params}
#' @description calculate density of parameters
#' @param x an object
#' @param ... additional arguments
#' @return density of parameters under x
#' @export
dparams <- function(x, ...)
    UseMethod("dparams", x)

#' @title Density of \code{params}
#' @description computing density of \code{params} under an \code{sbm} object
#' @param x an \code{sbm} object
#' @param mod a model list with a \code{params} element
#' @param ... additional arguments for \code{dparams.params}
#' @return density of parameters in \code{sbm} under mod
#' @export
dparams.sbm <- function(x, mod, ...)
    dparams.params(x$params, mod$params, ...)

#' @title Density of \code{params}
#' @description computing density of \code{params} under a \code{parammod}
#' @param x a \code{params} object
#' @param parammod a \code{parammod} object
#' @param log return log density?
#' @param ... additional parameters for \code{parammod$df} density function
#' @return density of parameters
#' @export
dparams.params <- function(x, parammod, log=TRUE, ...){
    out <- parammod$df(x$theta0, parammod$alpha, log=log,...) + sum(parammod$df(x$thetak, parammod$beta, log=log,...))
    if(!log)
        out <- exp(out)
    out
}

#' @title  plot \code{params}
#' @description barplot of parameter values in a \code{params} object
#' @param x a \code{params} object
#' @param ... additional graphical arguments
#' @export
plot.params <- function(x, ...){
    graphics::barplot(zapsmall(c(x$theta0, x$thetak)),...)
}

#' @title Printer
#' @description print parameters
#' @param x a \code{params} object
#' @param ... additional arguments for formatting (unused)
#' @export
print.params <- function(x,...)
    cat(format(x, ...), "\n")

format.params <- function(params,...)
    c("Params object\nkappa = ", params$kappa, "\ntheta0 =", params$theta0,"\nthetak:\n", params$thetak)

#' @title is \code{params} check
#' @description check if an R object is a \code{params} object
#' @param x an R object
#' @return TRUE if x is a \code{params} object
#' @export
is.params <- function(x)
    inherits(x, "params")

#' Make a matrix of parameters
#' @title Parameter Matrix maker
#' @export
#' @param x object for dispatch
#' @param ... additional arguments for method
parammat <- function(x, ...)
    UseMethod("parammat", x)

#' @title Parameter Matrix maker
#' @description Make a matrix of parameters from a matrix of block assignments
#' @param zleft block assignments on the left
#' @param zright block assignments on the right
#' @param params the parameters object
#' @param ... (unused)
#' @return a matrix of parameters of size |\code{zleft}| x |\code{zright}|
parammat.matrix <- function(zleft, zright, params, ...){
    p <- parammat(params, dim(zleft)[1])
    out <- array(0, c(params$dtheta, dim(zleft)[2], dim(zright)[2]))
    for(d in 1:params$dtheta)
        out[d,,] <- t(zleft) %*% p[d,,] %*% zright
    out
}

#' @title Parameter Matrix maker
#' @description Make a matrix of parameters from a \code{blocks} and \code{params} object
#' @param x a \code{blocks} object
#' @param params a \code{params} object
#' @param ... (unused)
#' @return an \code{NxN} matrix \code{P}, with \code{P[i,j] = } the parameter governing edge \code{ij}
#' @export
parammat.blocks <- function(x, params, ...)
    parammat(zmat(x), zmat(x), params)

#' @title Parameter Matrix maker
#' @description Make a matrix of parameters from a \code{params} object
#' @param x a \code{params} object
#' @param kappa - number of blocks to compute for - needed in computing the sequential assignment during split move
#' @param ... (unused)
#' @return a matrix of parameters
#' @export
parammat.params <- function(x, kappa, ...){
    ## possible to extend passed the number of parameters
    ## so we can look at empty blocks connecting with \theta_0
    if(missing(kappa))
        kappa <- x$kappa
    out <- array(0, c(x$dtheta, kappa, kappa))
    if(kappa > 1){
        for(d in 1:x$dtheta){
            out[d,,] <- matrix(x$theta0[d], kappa, kappa)
            diag(out[d,,])[1:x$kappa] <- x$thetak[,d, drop=F]
        }
    } else{
        out[,,] <- x$thetak
    }
    out
}

#' @title Parameter Matrix maker
#' @description Make a matrix of parameters from an \code{sbm} object
#' @param x an \code{sbm} object
#' @param ... (unused)
#' @return a matrix of parameters
#' @export
parammat.sbm <- function(x, ...)
    parammat(x$blocks, x$params)

#' @title Parameter model for Beta
#' @description beta parameter model
#' @param a0 \code{theta_0 ~ Beta(a0,a1)}
#' @param a1 \code{theta_0 ~ Beta(a0,a1)}
#' @param b0 \code{theta_k ~ Beta(b0,b1)}
#' @param b1 \code{theta_k ~ Beta(b0,b1)}
#' @return a \code{parammod}
#' @importFrom stats dbeta rbeta
#' @export
param_beta <- function(a0, a1, b0, b1){
    parammod(
        c(a0, a1),
        c(b0, b1),
        function(n, p, ...) stats::rbeta(n, p[1], p[2], ...),
        function(x, p, log, ...) stats::dbeta(x, p[1], p[2], log=log, ...),
        function(x) log(x) - log(1-x),
        function(x) 1/(1 + exp(-x)),
        function(x) -log(x)-log(1-x)
    )
}

#' gamma parameter model
#' @param a0 \code{theta_0 ~ Gamma(a0,a1)}
#' @param a1 \code{theta_0 ~ Gamma(a0,a1)}
#' @param b0 \code{theta_k ~ Gamma(b0,b1)}
#' @param b1 \code{theta_k ~ Gamma(b0,b1)}
#' @return a \code{parammod}
#' @importFrom stats rgamma dgamma
#' @export
param_gamma <- function(a0, a1, b0, b1){
    parammod(
        c(a0, a1),
        c(b0, b1),
        function(n, p, ...) stats::rgamma(n, p[1], p[2], ...),
        function(x, p, log, ...) stats::dgamma(x, p[1], p[2], log=log, ...),
        function(x) log(x),
        function(x) exp(x),
        function(x) -log(x)
    )
}

#' @title Parameter model for Normal Model
#' @description Normal parameter model:
#' \code{theta_0 = (mu0, sigma0)}
#' \code{theta_k = (muk, sigmak)}
#' @param a0,a1 \code{mu0 ~ Normal(a0,a1)}
#' @param b0,b1 \code{sig_0 ~ Gamma(b0,b1)}
#' @param c0,c1 \code{muk ~ Normal(c0, c1)}
#' @param d0,d1 \code{sig_k ~ Gamma(d0,d1)}
#' @return \code{parammod} representing Normal distributed parameters
#' @importFrom stats rgamma dgamma rnorm dnorm
#' @export
param_norm <- function(a0, a1, b0, b1, c0, c1, d0, d1){
    parammod(
        c(a0, a1, b0, b1),
        c(c0, c1, d0, d1),
        function(n, p, ...) cbind(stats::rnorm(n, p[1], p[2], ...), stats::rgamma(n, p[3], p[4], ...)),
        function(x, p, log, ...) sum(stats::dnorm(x[,1], p[1], p[2], log=log, ...), stats::dgamma(x[,2], p[3], p[4], log=log, ...)),
        function(x) cbind(x[1], log(x[2])),
        function(x) cbind(x[1],exp(x[2])),
        function(x) -log(x[2])
    )
}

#' @title Parameter model for Negative Binomial
#' @description Negative Binomial parameter model:
#' \code{theta_0 = (mu0, sigma0)}
#' \code{theta_k = (muk, sigmak)}
#' @param a0,a1 \code{mu0 ~ Gamma(a0,a1)}
#' @param b0,b1 \code{sig_0 ~ Beta(b0,b1)}
#' @param c0,c1 \code{muk ~ Gamma(c0, c1)}
#' @param d0,d1 \code{sig_k ~ Beta(d0,d1)}
#' @return \code{parammod} representing Negative-Binomial distributed parameters
#' @importFrom stats rgamma dgamma rbeta dbeta
#' @export
param_nbin <- function(a0, a1, b0, b1, c0, c1, d0, d1){
    parammod(
        c(a0, a1, b0, b1),
        c(c0, c1, d0, d1),
        function(n, p, ...) cbind(stats::rgamma(n, p[1], p[2], ...), stats::rbeta(n, p[3], p[4], ...)),
        function(x, p, log, ...) stats::dgamma(x[,1], p[1], p[2], log=log, ...) + stats::dbeta(x[,2], p[3], p[4], log=log, ...),
        function(x) cbind(log(x[1]), logit(x[2])),
        function(x) cbind(exp(x[1]), logistic(x[2])),
        function(x) -log(x[1])-log(x[2])-log(1-x[2])
    )
}
