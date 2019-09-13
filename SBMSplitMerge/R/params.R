#' Class parammod
#' @param alpha parameters for theta_0
#' @param beta parameters for theta_k
#' @param rf random function to draw parameters (the prior) rf(n, x) draws n parameters with hyperparameters x
#' @param df density function for parameters df(theta, x) density of theta with hyperparameters x
#' @param t map taking parameter space to real line
#' @param invt map taking real line to parameter space
#' @param loggradt log of the gradient of the map t
#' @return a parammod object
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

#' Class params
#' @param theta0 between block parameter
#' @param thetak within block parameters
#' @return a params object
#' @export
params <- function(theta0, thetak){
    thetak <- as.matrix(thetak)
    theta0 <- as.matrix(theta0)
    out <- list(theta0 = theta0, thetak = thetak,
                kappa = nrow(thetak), dtheta = ncol(thetak))
    class(out) <- c(class(out), "params")
    out
}

#' generate parameters
#' @param kappa number of blocks
#' @param PM a parammod object
#' @return a params object
#' @export
rparams <- function(kappa, PM){
    theta0 <- PM$rf(1, PM$alpha)
    thetak <- PM$rf(kappa, PM$beta)
    params(theta0, thetak)
}

#' calculate density of parameters
#' @param x an object
#' @return density of parameters under x
#' @export
dparams <- function(x, ...)
    UseMethod("dparams", x)

#' calculate density of parameters
#' @param sbm an sbm object
#' @param mod a model object
#' @return density of parameters in sbm under mod
#' @export
dparams.sbm <- function(sbm, mod, ...)
    dparams.params(sbm$params, mod$params, ...)

#' calculate density of parameters
#' @param params a params object
#' @param PM a parammod object
#' @param log return log density?
#' @param ... additional parameters for PM$df density funtcion
#' @return density of parameters in sbm under mod
#' @export
dparams.params <- function(Params, PM, log=TRUE, ...){
    out <- PM$df(Params$theta0, PM$alpha, log=log,...) + sum(PM$df(Params$thetak, PM$beta, log=log,...))
    if(!log)
        out <- exp(out)
    out
}

#'@export
plot.params <- function(params, ...){
    barplot(zapsmall(c(params$theta0, params$thetak)),...)
}

#' @export
print.params <- function(params,...)
    cat(format(params, ...), "\n")

format.params <- function(params,...)
    c("Params object\nkappa = ", params$kappa, "\ntheta0 =", params$theta0,"\nthetak:\n", params$thetak)

#' @export
is.params <- function(params)
    inherits(params, "params")

#' Make a matrix of parameters
parammat <- function(x, ...)
    UseMethod("parammat", x)

#' Helper function -  Make a matrix of parameters
#' @param zleft block assignments on the left
#' @param zright block assignments on the right
#' @param params the parameters object
#' @return a matrix of parameters of size |zleft| x |zright|
parammat.matrix <- function(zleft, zright, params, ...){
    p <- parammat(params, dim(zleft)[1])
    out <- array(0, c(params$dtheta, dim(zleft)[2], dim(zright)[2]))
    for(d in 1:params$dtheta)
        out[d,,] <- t(zleft) %*% p[d,,] %*% zright
    out
}

#' Helper function -  Make a matrix of parameters
#' @param blocks a blocks object
#' @param params the parameters object
#' @return a matrix of parameters of size NxN where parammet(i,j) is the parameter governing edge i,j under the block assignments in blocks
parammat.blocks <- function(blocks, params, ...)
    parammat(zmat(blocks), zmat(blocks), params)

#' Helper function -  Make a matrix of parameters
#' @param params the parameters object
#' @param kappa - number of blocks to compute for - needed in computing the sequential assignmetn during split move
#' @return a matrix of parameters
parammat.params <- function(Params, kappa){
    ## possible to extend passed the number of parameters
    ## so we can look at empty blocks connecting with \theta_0
    if(missing(kappa))
        kappa <- Params$kappa
    out <- array(0, c(Params$dtheta, kappa, kappa))
    if(kappa > 1){
        for(d in 1:Params$dtheta){
            out[d,,] <- matrix(Params$theta0[d], kappa, kappa)
            diag(out[d,,])[1:Params$kappa] <- Params$thetak[,d, drop=F]
        }
    } else{
        out[,,] <- Params$thetak
    }
    out
}

#' Helper function -  Make a matrix of parameters
#' @param sbm an SBM object
parammat.sbm <- function(sbm)
    parammat(sbm$blocks, sbm$params)

#' beta parameter model
#' @param a0 theta_0 ~ Beta(a0,a1)
#' @param a1 theta_0 ~ Beta(a0,a1)
#' @param b0 theta_k ~ Beta(b0,b1)
#' @param b1 theta_k ~ Beta(b0,b1)
#' @return parammod representing beta distributed parameters
#' @export
param_beta <- function(a0, a1, b0, b1){
    parammod(
        c(a0, a1),
        c(b0, b1),
        function(n, p, ...) rbeta(n, p[1], p[2], ...),
        function(x, p, log, ...) dbeta(x, p[1], p[2], log=log, ...),
        function(x) log(x) - log(1-x),
        function(x) 1/(1 + exp(-x)),
        function(x) -log(x)-log(1-x)
    )
}

#' gamma parameter model
#' @param a0 theta_0 ~ Gamma(a0,a1)
#' @param a1 theta_0 ~ Gamma(a0,a1)
#' @param b0 theta_k ~ Gamma(b0,b1)
#' @param b1 theta_k ~ Gamma(b0,b1)
#' @return parammod representing gamma distributed parameters
#' @export
param_gamma <- function(a0, a1, b0, b1){
    parammod(
        c(a0, a1),
        c(b0, b1),
        function(n, p, ...) rgamma(n, p[1], p[2], ...),
        function(x, p, log, ...) dgamma(x, p[1], p[2], log=log, ...),
        function(x) log(x),
        function(x) exp(x),
        function(x) -log(x)
    )
}

#' Normal parameter model
#' theta_0 = (mu0, sigma0)
#' theta_k = (muk, sigmak)
#' @param a0,a1 mu0 ~ Normal(a0,a1)
#' @param b0,b1 sig_0 ~ Gamma(b0,b1)
#' @param c0,c1 muk ~ Normal(c0, c1)
#' @param d0,d1 sig_k ~ Gamma(d0,d1)
#' @return parammod representing Normal distributed parameters
#' @export
param_norm <- function(a0, a1, b0, b1, c0, c1, d0, d1){
    parammod(
        c(a0, a1, b0, b1),
        c(c0, c1, d0, d1),
        function(n, p, ...) cbind(rnorm(n, p[1], p[2], ...), rgamma(n, p[3], p[4], ...)),
        function(x, p, log, ...) sum(dnorm(x[,1], p[1], p[2], log=log, ...), dgamma(x[,2], p[3], p[4], log=log, ...)),
        function(x) cbind(x[1], log(x[2])),
        function(x) cbind(x[1],exp(x[2])),
        function(x) -log(x[2])
    )
}

#' Negative Binomial parameter model
#' theta_0 = (mu0, sigma0)
#' theta_k = (muk, sigmak)
#' @param a0,a1 mu0 ~ Gamma(a0,a1)
#' @param b0,b1 sig_0 ~ Beta(b0,b1)
#' @param c0,c1 muk ~ Gamma(c0, c1)
#' @param d0,d1 sig_k ~ Beta(d0,d1)
#' @return parammod representing NegativeBinomial distributed parameters
#' @export
param_nbin <- function(a0, a1, b0, b1, c0, c1, d0, d1){
    parammod(
        c(a0, a1, b0, b1),
        c(c0, c1, d0, d1),
        function(n, p, ...) cbind(rgamma(n, p[1], p[2], ...), rbeta(n, p[3], p[4], ...)),
        function(x, p, log, ...) dgamma(x[,1], p[1], p[2], log=log, ...) + dbeta(x[,2], p[3], p[4], log=log, ...),
        function(x) cbind(log(x[1]), logit(x[2])),
        function(x) cbind(exp(x[1]), logistic(x[2])),
        function(x) -log(x[1])-log(x[2])-log(1-x[2])
    )
}
