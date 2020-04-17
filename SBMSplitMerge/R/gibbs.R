#' Gibbs step for drawing new block assignments
#' @param SBM an \code{sbm} object
#' @param Edges an \code{edges} object
#' @param Mod a model list
#' @return updated \code{sbm} object
drawblocks.gibbs <- function(SBM, Edges, Mod){
    for(i in 1:SBM$numnodes)
        SBM <- drawblock.gibbs(i, SBM, Edges, Mod)
    SBM
}

#' Gibbs step for drawing new block assignments
#' @param i the node to reassign
#' @param SBM an \code{sbm} object
#' @param Edges an \code{edges} object
#' @param Mod a model list
#' @return updated \code{sbm} object
drawblock.gibbs <- function(i, SBM, Edges, Mod){
    p <- condprior(SBM$blocks, Mod$blocks, i) + nodelike(Edges, SBM, Mod, i)[1:SBM$kappa]     ## assume that we can't make a new block here for rjmcmc and vanilla gibbs
    updateblock(SBM, rcat(1, normaliselogs(p)), i)
}
