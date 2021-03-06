# SBM SplitMerge
This repo contains the source for the SBMSplitMerge R package for performing inference on the number of blocks and the edge-state parameters in a Generalised Stochastic Block Model using a split-merge sampler.
This allows for inference in non-conjugate edge-state models.

This repo also contains the examples and data application scripts to reproduce the plots in the accompanying paper on [ArXiV https://arxiv.org/abs/1909.09421](https://arxiv.org/abs/1909.09421).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3407683.svg)](https://doi.org/10.5281/zenodo.3407683)


To install from github using ``R`` issue:
```r
devtools::install_github("ludkinm/SBMSplitMerge", subdir="SBMSplitMerge", build_vignettes = TRUE)
```
After installation, see the vignette for example usage:
```R
vignette("Weibull-edges")
```
