# L4AssessRPackage
 State space biomass dynamics models incorporating environmental effects

"L4Assess" contains a series of 'wrapper' R functions for applying dynamic, state space biomass dynamics models, implemented in template model builder (TMB). The various routines
that have been implemtned include various modelling options, including specified production functions (Schaefer, Fox or Pella-Tomlinson), use of  environmental index data 
to account for environmental effects on stock dynamics and/or catchability, allow for depensatory stock dynamics, and if available, fit to estimates of absolute abundance 
from fishery-independent survey data. Stock dynamics are modelled using the Schaefer, Fox or Pella-Tomlinson production functions. These models are typically fitted to 
time series data for annual catches (total stock removals) and one or more indices of stock abundance, such as fishery-independent or standardised, fishery dependent 
catch per unit effort (CPUE) data. 

Note that biomass dynamics models implictly assume, among other things, that fish become mature at the same time that they become selected into the fishery coincide (and
also that this remains constant over time). Violation of these assumptions can affect results. If sufficient data are available, core sophisticated dynamic  approaches 
(i.e. integrated models) which model these different processes explictly are recommended, particularly if this assumption of biomass dynamics models is unlikely to be met.

To run the available TMB biomass dynamics models using the L4Assess R package, currently, it is necessary to have the TMB R package installed, and to download 
the TMB code (SSBDM.cpp) for the TMB biomass dynamics models from github ("SAlexHesp/L4AssessRPackage"). It may also be helpful to download the example R script 
("SSBDM code to run model.R") and example data ("SSBDM CS crab data.csv") from the same github site.

As with other packages in github, to install L4Assess, first ensure you have the 'devtools' package, otherwise this can be installed using install.packages("devtools').
Then, use: 

library(devtools)

devtools::install_github("SAlexHesp/L4AssessRPackage", build_vignettes=TRUE)

If you already have a version of the L4Assess package installed but wish to update, I suggest using the line of code below rather than the one above 
to ensure the updated version is installed.
devtools::install_github("SAlexHesp/L4AssessRPackage", build_vignettes=TRUE, force=TRUE)
