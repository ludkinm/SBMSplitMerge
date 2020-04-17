## Misc
normaliselogs <- function(logp){
    a <- max(logp)
    exp(logp - a - log(sum(exp(logp - a))))
}

logit <- function(p)
    log(p)-log(1-p)

logistic <- function(x)
    1/(1+exp(-x))

xlogx <- function(x)
    ifelse(x==0, x, x*log(x))

#' Draw n random variables from Dricihlet(gam) - gam a vector of length K
#' @param n number of variates to draw
#' @param gam a vector of concentration parameters
#' @return x = Matrix(n,K)
#' @importFrom stats rgamma
#' @export
rdirichlet <- function (n, gam) {
    l <- length(gam)
    x <- matrix(stats::rgamma(l * n, gam), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    x/as.vector(sm)
}

#' Density of X~dirichlet(gam)
#' @param x random variable in the d-dimensional simplex
#' @param gam a length d probability vector
#' @param log return the log-probability instead?
#' @return the density
#' @export
ddirichlet <- function (x, gam, log=FALSE) {
       if(length(gam) != (nc <- ncol(x)))
        gam <- rep(gam,nc)
    p <- sum(lgamma(gam)) - lgamma(sum(gam)) + rowSums(gam * log(x))
    if(!log)
        p <- exp(p)
    p
}

#' Random draw from Categorical(p)
#' @param n number of draws
#' @param p a length-d probability vector
#' @param replace should the categories be replaced? If so n < p required
#' @return a length-d vector x, where \eqn{x[i] \in \{1,...,d\}}
#' @export
rcat <- function(n, p, replace=TRUE)
    sample(1:length(p), n, replace, p)

#' V-measure
#' @param z input vector
#' @param truez reference vector
#' @param beta vmeasure parameter beta=1 gives equal weight to homogeneity and completeness
#' @return v-measure of z against truez
#' @export
vmeasure <- function(z, truez, beta=1){
    C <- z
    B <- truez
    A <- table(B,C)
    N <- sum(A)

    HB  <- log(N) - sum(xlogx(rowSums(A)))/N
    HC  <- log(N) - sum(xlogx(colSums(A)))/N
    HBC <- sum(xlogx(colSums(A)))/N - sum(xlogx(A))/N
    HCB <- sum(xlogx(rowSums(A)))/N - sum(xlogx(A))/N

    h <- 1 - HBC/HB
    c <- 1 - HCB/HC
    v <- (beta+1)*h*c/(beta*h+c)
    v
}

#' Adjusted Rand Index
#' @param z input vector
#' @param truez reference vector
#' @return ARI of z against truez
#' @export
ARI <- function(z, truez){
    N <- length(z)
    n <- table(z, truez)
    a <- sum(choose(colSums(n),2))
    b <- sum(choose(rowSums(n),2))
    eind <- a*b/choose(N,2)
    top <- sum(choose(n,2)) - eind
    bot <- 0.5*(a+b) - eind
    top/bot
}


#' gets command line args to numerics
#' @export
#' @param cmdargs a vector of strings probably from a call to commandArgs
#' @param i the index to extract
#' @param default a default value to use if not in cmdargs
get_arg <- function(cmdargs, i, default){
    tmp <- as.numeric(cmdargs[i])
    if(!is.na(tmp)){
        return(tmp)
    } else{
        return(default)
    }
}
