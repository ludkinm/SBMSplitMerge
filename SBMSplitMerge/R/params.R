#' @title \code{params} S3 object
#' @description make a \code{params} object from the between-block parameter \code{theta0} and a vector of within block parameters \code{thetak}
#' @param theta0 between block parameters - a vector of length `dimension of theta`
#' @param thetak within block parameters - a matrix with ncol=kappa and nrow=dimension of theta
#' @return a \code{params} object
#' @export
#' @examples
#' p <- params(0.1, c(0.2,0.4,0.5))
#' p
params <- function(theta0, thetak){
    theta0 <- c(theta0)
    thetak <- as.matrix(thetak)
    out <- list(
        theta0 = theta0,
        thetak = thetak,
        dimtheta = ncol(thetak),
        kappa = nrow(thetak)
    )
    class(out) <- c(class(out), "params")
    out
}

print.params <- function(x,...){
    cat("Params object\nkappa: ", x$kappa, "\ntheta0:")
    print(x$theta0)
    cat("\nthetak:\n")
    print(x$thetak)
}

is.params <- function(x)
    inherits(x, "params")

#' @importFrom graphics barplot
plot.params <- function(x, xlab="Block", ylab="Value", beside=TRUE, ...){
    graphics::barplot(t(rbind(x$theta0, x$thetak)), beside=beside, xlab=xlab, ylab=ylab, ...)
}
