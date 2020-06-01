Adressed issues raised by Jelena Saf:
1) Shortened the title to a maximum of 65 characters
Specifically replaced Stochastic Block Model with SB and defined in this in `Description`.

2) Added these missing Rd-tags:
	  blockmat.Rd: \value
	  plot.blocks.Rd: \value
	  plot.sbm.Rd: \value
	  plotpostpairs.Rd: \value
	  splitparams.Rd: \value
	  updateblock.Rd: \value

Passed checks:
- R CMD check
- devtools::check(cran=TRUE)
- devtools::check_rhub():
- - https://builder.r-hub.io/status/SBMSplitMerge_1.1.0.tar.gz-c8a9ef32a13d428ba82b1f8b7cd9d22a
- - https://builder.r-hub.io/status/SBMSplitMerge_1.1.0.tar.gz-e7a4b7c1cc6a421698ce99044148ffc6
- - https://builder.r-hub.io/status/SBMSplitMerge_1.1.0.tar.gz-25affe58944348a089588bcc44304d3f
-  devtools::check_win_devel()
- - https://win-builder.r-project.org/9fPZ1lBCRWta/
