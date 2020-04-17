#' create a blocks object
#' @param z vector of block labels
#' @param kappa maximum number of blocks
#' @return a blocks object
#' @examples
#' b <- blocks(c(1,1,2,3,4,4), 4)
#' @export
blocks <- function(z, kappa){
    z <- as.factor(z)
    if(missing(kappa))
        kappa <- nlevels(z)
    levels(z) <- 1:kappa
    out <- list(
        z = z,
        kappa = kappa,
        sizes = table(z),
        numnodes = length(z)
    )
    class(out) <- c(class(out), "blocks")
    out
}

#' create a block-model object
#' @param gamma parameters for the block model
#' @param fixkappa Logical - is kappa fixed or can it vary
#' @param rblocks function(n, gamma) samples node memberships from the block model, takes gamma as a parameter
#' @param logdblocks log density for the number of blocks distribution
#' @param dcond conditional density for the block assignments
#' @return a \code{blockmod} object
#' @export
blockmod <- function(gamma, fixkappa, rblocks, logdblocks, dcond){
    out <- list(
        gamma = gamma,
        fixkappa = fixkappa,
        rblocks = rblocks,
        logdblocks = logdblocks,
        dcond = dcond
    )
    class(out) <- c(class(out), "blockmod")
    out
}

#' plots a block object
#' @param x a blocks object to plot
#' @param col colours for the plot
#' @param xaxt override \code{image} parameters
#' @param yaxt override \code{image} parameters
#' @param xlab override \code{image} parameters
#' @param ylab override \code{image} parameters
#' @param ... additional parameters for \code{image}
#' @importFrom graphics image plot
#' @importFrom grDevices rainbow
#' @export
plot.blocks <- function(x, col, xaxt='n', yaxt='n', xlab="", ylab="", ...){
    tmp <- zmat(x)
    if(missing(col)){
        if(x$kappa > 1){
            col <- c(0, grDevices::rainbow(x$kappa))
        } else{
            col <- grDevices::rainbow(x$kappa)
        }
    }
    n <- 1:x$numnodes-0.5
    graphics::image(n, n, t(tmp) %*% ((1:x$kappa) * tmp), col=col,
          xaxt=xaxt, yaxt=yaxt, xlab=xlab, ylab=ylab, ...)
}

#' @rdname plot.blocks
#' @export
image.blocks <- plot.blocks

#' print a \code{blocks} object
#' @param x a \code{blocks} object
#' @param ... additional arguments for formatting
#' @export
print.blocks <- function(x, ...)
    cat(format(x, ...), "\n")

format.blocks <- function(blocks, ...)
    c("blocks object\nkappa = ", blocks$kappa, "\nNumber of nodes =", blocks$numnodes,"\nblock sizes:\n", blocks$sizes)

#' check if an \code{R} object is a \code{blocks} object
#' @param x an \code{R} object to dispatch on
#' @param ... additional arguments
#' @export
is.blocks <- function(x, ...)
    inherits(x, "blocks")

#' converts x to a matrix of block assignments
#' @param x object for dispatch
#' @param ... additional arguments for method
zmat <- function(x, ...)
    UseMethod("zmat", x)

#' converts x to a matrix of block assignments
#' @param x a vector of node-to-block assignments
#' @param K number of blocks
#' @return matrix with K rows and a 1 at (k,i) if node i is in block k under x
zmat.numeric <- function(x, K)
    diag(K)[ , x, drop=FALSE]

#' converts x to a matrix of block assignments
#' @param x a factor of node-to-block assignments
#' @param K number of blocks
#' @return matrix with K rows and a 1 at (k,i) if node i is in block k under x
zmat.factor <- zmat.numeric

#' converts blocks to a matrix of block assignments
#' @param blocks a blocks object
#' @param kappa number of blocks
#' @return matrix with K rows and a 1 at (k,i) if node i is in block k under the blocks object
zmat.blocks <- function(blocks, kappa){
    if(missing(kappa))
        kappa <- blocks$kappa
    zmat(blocks$z, kappa)
}

#' converts the block assignment under \code{SBM} to a matrix of block assignments
#' @param SBM an \code{sbm} object
#' @return matrix with K rows and a 1 at (k,i) if node i is in block k under the \code{sbm} model
zmat.sbm <- function(SBM)
    zmat(SBM$blocks)

#' density for the block assignments in x
#' @param x object for dispatch
#' @param ... additional arguments for method
#' @export
dblocks <- function(x, ...)
    UseMethod("dblocks", x)

#' density for the block assignments
#' @param x an \code{sbm} model which contains the block state
#' @param mod a list containing an element \code{mod$blocks} a \code{blockmod}
#' @param ... additional arguments for \code{dblocks.blocks}
#' @return The density for the block assignments in \code{sbm} under \code{mod}
#' @export
dblocks.sbm <- function(x, mod, ...)
    dblocks.blocks(x$blocks, mod$blocks, ...)

#' density for the block assignments
#' @param x a \code{blocks} object containing block assignments
#' @param blockmod a \code{blockmod} object containing a model for the block assignments
#' @param ... additional arguments for \code{blockmod$logdblocks}
#' @return The density for the block assignments in \code{x} under \code{blockmod}
#' @export
dblocks.blocks <- function(x, blockmod, ...)
    blockmod$logdblocks(x, blockmod$gamma, ...)

#' random block assignment generation
#' @param n number of nodes
#' @param blockmod a \code{blockmod} object
#' @param sorted sort the nodes by block assignment for nice plots?
#' @param ... additional parameters for the \code{blockmod$rblocks} function
#' @return a blocks object
#' @export
rblocks <- function(n, blockmod, sorted=FALSE, ...){
    z <- blockmod$rblocks(n, blockmod$gamma, ...)
    if(sorted)
        z <- sort(z)
    blocks(z)
}

#' @title Conditional Prior
#' @description Calculate the conditional density for block assignment of node i under the model in \code{blockmod} and block assignments in blocks
#' @param blocks \code{blocks} object
#' @param blockmod a \code{blockmod} object
#' @param i a node index
#' @param ... additional arguments for \code{blockmodel$dcond}
#' @return conditional density for block assignment of node i under the model in \code{blockmod} and block assignments in blocks
#' @export
condprior <- function(blocks, blockmod, i, ...)
    blockmod$dcond(blocks, blockmod$gamma, i, ...)


#' Update the block assignment of a node
#' @param x object for dispatch
#' @param ... additional arguments for method
updateblock <- function(x,...)
    UseMethod("updateblock")

#' Update the block assignment of a node
#' @param SBM an \code{sbm} object
#' @param newblock the new block for node i
#' @param i the node to update
#' @return new \code{sbm} object
updateblock.sbm <- function(SBM, newblock, i){
    newblocks <- updateblock(SBM$blocks, newblock,i)
    sbm(newblocks, SBM$params)
}

#' Update the block assignment of a node
#' @param blocks a \code{blocks} object
#' @param newblock the new block for node i
#' @param i the node to update
#' @return new \code{blocks} object
updateblock.blocks <- function(blocks, newblock, i){
    tmp <- blocks$z
    tmp[i] <- newblock
    blocks(tmp)
}

#' Multinomial block assignment
#' @param gamma parameter for Dirichlet component \code{Dir(gam, gam, .., gam)}
#' @return a block model representing a \code{Multinomial(gamma)} distribution
#' @export
multinom <- function(gamma){
    blockmod(
        gamma
       ,
        TRUE
       ,
        function(n, gamma, kappa){
            omega <- rdirichlet(1, rep(gamma[1], kappa))
            rcat(n, c(omega))
        }
       ,
        function(blocks, gam){
            gam <- gam[1]
            kap <- blocks$kappa
            out <- lgamma(kap * gam) - kap * lgamma(gam) + sum(lgamma(blocks$sizes + gam)) - lgamma(blocks$numnodes + kap*gam)
            out
        }
       ,
        function(Blocks, gamma, i){
            ## log(number of nodes per block without i + gamma)
            Nik <- Blocks$sizes - tabulate(Blocks$z[i], Blocks$kappa)
            log(Nik + gamma[1]) - log(sum(Nik + Blocks$kappa*gamma[1]))
        }
    )
}

#' Dirichlet Multinomial Allocation model
#' @param gamma parameter for Dirichlet component \code{Dir(gam, gam, ...)}
#' @param delta parameter for \code{Poison(delta)} - distribution on number of blocks
#' @return a block model representing a \code{dma(gamma, delta)} distribution
#' @export
dma <- function(gamma, delta){
    blockmod(
        c(gamma, delta),
        FALSE,
        function(n, gamma, ...){
            kappa <- stats::rpois(1, gamma[2]) + 1
            omega <- rdirichlet(1, rep(gamma[1], kappa))
            rcat(n, c(omega))
        },
        function(blocks, gamma){
            gam <- gamma[1]
            del <- gamma[2]
            kap <- blocks$kappa
            stats::dpois(kap-1, del, log=T) + lgamma(kap * gam) - kap * lgamma(gam) + sum(lgamma(blocks$sizes + gam)) - lgamma(blocks$numnodes + kap*gam)
        },
        function(Blocks, gamma, i){
            ## For a fixed number of blocks, the conditional distribution that node
            ## i belongs to block k given gamma (log scale)
            kappa <- Blocks$kappa
            Mk <- Blocks$sizes - c(zmat(Blocks$z[i], kappa))
            return( log(Mk + gamma[1]) - log(sum(Mk) + kappa*gamma[1]) )
        }
    )
}

#' compute log gamma with \code{lgamma(0) = -Inf}
#' @param x numeric value to compute log gamma
mylgamma <- function(x)
    sapply(x, function(x) ifelse(x==0, -Inf, lgamma(x)))


#' Chinese Restaurant Process model
#' @param gamma concentration parameter
#' @return a block model representing a \code{CRP(gamma)} distribution
#' @export
crp <- function(gamma){
    blockmod(
        gamma
       ,
        FALSE
       ,
        function(n, gamma, ...){
            z <- c(1, rep(0,n-1))
            blocksizes <- 1
            if(n > 1)
                for(i in 2:n){
                    probs <- c(tabulate(z), gamma)
                    z[i] <- sample(1:length(probs), 1, FALSE, probs)
                }
            z
        },
        function(blocks, gamma){
            gam <- gamma[1]
            kap <- blocks$kappa
            lgamma(gam) + kap * gam  - lgamma(blocks$numnodes + gam) + sum(mylgamma(blocks$sizes))
        },
        function(blocks, gamma, i)
            log(c(tabulate(blocks$z[-i], blocks$kappa), gamma))
    )
}
