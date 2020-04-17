#' Gibbs step for drawing new block assignments
#' @param Sbm an sbm object
#' @param Edges an edges object
#' @param Mod a model object
#' @return updated sbm object
drawblocks.gibbs <- function(Sbm, Edges, Mod){
    for(i in 1:Sbm$numnodes)
        Sbm <- drawblock.gibbs(i, Sbm, Edges, Mod)
    Sbm
}

#' Gibbs step for drawing new block assignments
#' @param i the node to reassign
#' @param Sbm an sbm object
#' @param Edges an edges object
#' @param Mod a model object
#' @return updated sbm object
drawblock.gibbs <- function(i, Sbm, Edges, Mod){
    p <- condprior(Sbm$blocks, Mod$blocks, i) + nodelike(Edges, Sbm, Mod, i)[1:Sbm$kappa]     ## assume that we can't make a new block here for rjmcmc and vanilla gibbs
    updateblock(Sbm, rcat(1, normaliselogs(p)), i)
}
