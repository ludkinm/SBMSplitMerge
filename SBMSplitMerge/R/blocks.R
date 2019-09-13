#' creat a blocks object
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
#' @param fixkappa bool - is kappa fixed or can it vary
#' @param rblocks function(n, gamma) samples node memberships from the block model, takes gamma as a parameter
#' @param logdblocks log density for the number of blocks distribution
#' @param dmarg marginal density for the block assignments
#' @return a bockmod object
#' @export
blockmod <- function(gamma, fixkappa, rblocks, logdblocks, dmarg){
    out <- list(
        gamma = gamma,
        fixkappa = fixkappa,
        rblocks = rblocks,
        logdblocks = logdblocks,
        dmarg = dmarg
    )
    class(out) <- c(class(out), "blockmod")
    out
}

#' plots a block object
#' @param blocks a blocks object to plot
#' @param col colors for the plot
#' @param xaxt='n' override 'image' parameters
#' @param yaxt='n' override 'image' parameters
#' @param xlab="" override 'image' parameters
#' @param ylab="" override 'image' parameters
#' @param ... additional parameters for 'image'
#' @export
plot.blocks <- image.blocks <- function(blocks, col, xaxt='n', yaxt='n', xlab="", ylab="", ...){
    tmp <- zmat(blocks)
    if(missing(col)){
        if(blocks$kappa > 1){
            col <- c(0, rainbow(blocks$kappa))
        } else{
            col <- rainbow(blocks$kappa)
        }
    }
    n <- 1:blocks$numnodes-0.5
    image(n, n, t(tmp) %*% ((1:blocks$kappa) * tmp), col=col,
          xaxt=xaxt, yaxt=yaxt, xlab=xlab, ylab=ylab, ...)
}

#' @export
print.blocks <- function(blocks,...)
    cat(format(blocks, ...), "\n")

#' @export
format.blocks <- function(blocks,...)
    c("blocks object\nkappa = ", blocks$kappa, "\nNumber of nodes =", blocks$numnodes,"\nblock sizes:\n", blocks$sizes)

#' @export
is.blocks <- function(blocks)
    inherits(blocks, "blocks")

#' converts x to a matrix of block assignments
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

#' converts the block assignment under sbm to a matrix of block assignments
#' @param sbm an sbm object
#' @return matrix with K rows and a 1 at (k,i) if node i is in block k under the sbm model
zmat.sbm <- function(sbm)
    zmat(sbm$blocks)

#' density for the block assignments in x
#' @export
dblocks <- function(x,...)
    UseMethod("dblocks", x)

#' density for the block assignments
#' @param sbm an sbm model which contains the block state
#' @param mod a model for the block assignments
#' @return The density for the block assignments in sbm under mod
#' @export
dblocks.sbm <- function(sbm, mod)
    dblocks.blocks(sbm$blocks, mod$blocks)

#' density for the block assignments
#' @param blocks a blocks object containing block assignments
#' @param blockmod a blockmod object containing a model for the block assignments
#' @return The density for the block assignments in blocks under blockmod
#' @export
dblocks.blocks <- function(blocks, blockmod, ...)
    blockmod$logdblocks(blocks, blockmod$gamma)

#' random block assignment generation
#' @param n number of nodes
#' @param BM a blockmod object
#' @param sorted sort the nodes by block assignment for nice plots?
#' @param ... additional parameters for the BM$rblocks function
#' @return a blocks object
#' @export
rblocks <- function(n, BM, sorted=FALSE, ...){
    z <- BM$rblocks(n, BM$gamma, ...)
    if(sorted)
        z <- sort(z)
    blocks(z)
}

#' Calculate the marginal density for block assignment of node i under the model in BM and block assignments in blocks
#' @param blocks blocks object
#' @param BM a blockmod object
#' @param i a node index
#' @return marginal density for block assignment of node i under the model in BM and block assignments in blocks
#' @export
margprior <- function(blocks, BM, i, ...)
    BM$dmarg(blocks, BM$gamma, i)


#' Update the block assignment of a node
updateblock <- function(x,...)
    UseMethod("updateblock")

#' Update the block assignment of a node
#' @param sbm an sbm object
#' @param newblock the new block for node i
#' @param i the node to update
#' @return new sbm object
updateblock.sbm <- function(sbm, newblock, i){
    newblocks <- updateblock(sbm$blocks, newblock,i)
    sbm(newblocks, sbm$params)
}

#' Update the block assignment of a node
#' @param blocks a blocks object
#' @param newblock the new block for node i
#' @param i the node to update
#' @return new blocks object
updateblock.blocks <- function(blocks, newblock, i){
    tmp <- blocks$z
    tmp[i] <- newblock
    blocks(tmp)
}

#' Multionoial block assignment
#' @param gamma parameter for dirichlet component Dir(gam, gam, .., gam)
#' @return a block model representing a multinomial(gamma) distribution
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

#' Dirichlet Multionoial Allocation model
#' @param gamma parameter for dirichlet component Dir(gam, gam, ...)
#' @param delta parameter for Poison(delta) - distribution on number of blocks
#' @return a block model representing a dma(gamma, delta) distribution
#' @export
dma <- function(gamma, delta){
    blockmod(
        c(gamma, delta),
        FALSE,
        function(n, gamma, ...){
            kappa <- rpois(1, gamma[2]) + 1
            omega <- rdirichlet(1, rep(gamma[1], kappa))
            rcat(n, c(omega))
        },
        function(blocks, gamma){
            gam <- gamma[1]
            del <- gamma[2]
            kap <- blocks$kappa
            dpois(kap-1, del, log=T) + lgamma(kap * gam) - kap * lgamma(gam) + sum(lgamma(blocks$sizes + gam)) - lgamma(blocks$numnodes + kap*gam)
        },
        function(Blocks, gamma, i){
            ## For a fixed number of blocks, the conditional distribution that node
            ## i belongs to block k given gamma (log scale)
            kappa <- Blocks$kappa
            Mk <- Blocks$sizes - c(zmat(Blocks$z[i], kappa))
            return( log(Mk + gamma[1]) - log(sum(Mk) + kappa*gamma[1]) )
            ## p   <- log( eta/(eta+gamma[2]) * c((Mk + gamma[1])/(M+eta*gamma[1]),0) + gamma[2]/(eta+gamma[2]) * beta(gamma[1] + 1, M+eta*gamma[1])/beta(gamma[1], eta*gamma[1])*c(1+Mk/gamma[1], 1))
            ## p
        }
    )
}

#' compute lgamma with lgamma(0) = -Inf
mylgamma <- function(x)
    sapply(x, function(x) ifelse(x==0, -Inf, lgamma(x)))


#' Chinese Restaurant Process model
#' @param gamma parameter for CRP(gam)
#' @return a block model representing a CRP(gamma) distribution
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
