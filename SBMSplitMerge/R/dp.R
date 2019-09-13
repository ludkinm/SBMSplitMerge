#' draw blocks in Dirichlet process sampler
#' @param currsbm current SBM state
#' @param edges an edges object
#' @param mod a model object
#' @return updated sbm
drawblocks.dp <- function(currsbm, Edges, Mod){
    for(i in 1:currsbm$numnodes)
        currsbm <- drawblock.dp(i, currsbm, Edges, Mod)
    currsbm
}

#' draw blocks in Dirichlet process sampler
#' @param i node to update
#' @param currsbm current SBM state
#' @param edges an edges object
#' @param mod a model object
#' @return updated sbm
drawblock.dp <- function(i, currsbm, Edges, Mod){
    ## current values
    propz <- as.numeric(as.character(currsbm$blocks$z))
    currb <- propz[i]
    pthetak <- currsbm$params$thetak
    ptheta0 <- currsbm$params$theta0
    pblocks <- currsbm$blocks

    ## probability calculations
    p0 <- nodelike(Edges, currsbm, Mod, i)
    p1 <- margprior(currsbm$blocks, Mod$blocks, i)
    p <- p0+p1

    ## choose new block
    propb <- rcat(1, normaliselogs(p))
    propz[i] <- propb

    ## Cases
    ## 1) (propb == currb) - do nothing i is assigned to its current block
    ## 2) (propb != currb) and size[currb] == 1 and propb == currkappa+1 - currb drops out of model propb added to model -> but labels dont matter so just simulate a new parameter
    ## 3) (propb != currb) and size[currb] > 1  and propb == currkappa+1 - i moves from currb to a new block - simulate a new parameter for the new block
    ## 4) (propb != currb) and size[currb] == 1 and propb != currkappa+1 - currb drops out of model
    ## 5) (propb != currb) and size[currb] > 1  and propb != currkappa+1 - i moves from currb to propb

    if( propb != currb ){
        ## then something changes...
        if( (currsbm$blocks$sizes[currb] == 1) && (propb == (currsbm$blocks$kappa+1)) ) {
            ## edge case - currb is a singleton propb is a new block:
            ## just resample parameter since labels dont matter
            pthetak[currb,] <- rparams(1, Mod$params)$thetak
        } else{
            ## otherwise we will move i
            pblocks <- blocks(propz)
            if( (currsbm$blocks$sizes[currb] == 1) && (propb != (currsbm$blocks$kappa+1)) ) {
                ## moving i leaves currb empty and will be removed from the sbm-state:
                pthetak <- pthetak[-currb,,drop=FALSE]
            } else if( (currsbm$blocks$sizes[currb] > 1) && (propb == (currsbm$blocks$kappa+1)) ) {
                ## Moving to a new block which needs a parameter
                pthetak <- rbind(pthetak, rparams(1, Mod$params)$thetak)
            }
        }
    }

    currsbm <- sbm(pblocks, params(ptheta0, pthetak))
    currsbm
}
