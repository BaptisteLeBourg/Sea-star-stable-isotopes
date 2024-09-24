# Sea-star-stable-isotopes

This repository provides the script used in the paper for "How do environmental conditions affect sea star trophic ecology? A circum-Antarctic assessment" by Le Bourg et al. (2024).

It contains two tables: 

IsotopeData.csv provides stable isotope data in sea stars from the Southern Ocean. The columns are:
- ExpeditionID: name of the sampling expedition
- StationID: name of the sampling station
- ScientificName: species name of the sampled sea star
- Date: date of the sampling
- Latitude: latitude of the sampling
- Longitude: longitude of the sampling
- region: name of the region where the sampling occurred
- Depth: depth at which the sea star was collected
- d13C: stable isotope values of carbon
- d15N: stable isotope values of nitrogen
- Trophic_group: trophic group of the collected sea star (?: unknown trophic group, Carnivore1: predator of active prey, Deposit-feeder: sediment feeder, Omnivore: omnivore, Pelagic-Carnivore1: predator of sessile prey, Pelagic-Carnivore2: predator of encrusting prey, Pelagic-Omnivore: pelagos-base omnivore, Pelagic-Suspension-feeder: suspension feeder)

EnvironmentalData.csv provides environmental data for each station where sea stars were sampled. The columns are:
- ExpeditionID: name of the sampling expedition
- StationID: name of the sampling station
- Date: date of the sampling
- Latitude: latitude of the sampling
- Longitude: longitude of the sampling
- seaice_prev_month: sea ice concentration during the month preceding the sampling
- days_since_melt: numbers of days between sea ice melting and the sampling. Values equal to 32765 corresponds to areas that were never recorded to be covered by sea ice.
- seaice_last_730: number of days during which sea ice concentration was higher than 85 % over a period of 730 days preceding the sampling
- chl_prev_month: surface chlorophyll concentration from the month preceding the sampling

Please, note that these are not the complete datasets used in the paper. Indeed, some of the data used in the paper were from previous studies and were thus not shared here (with the exception of the data from Michel et al., 2019, Scientific Reports 9: 8062). Consequently, the results displayed by the scripts will slightly differ from the ones shown in the paper.

Most data provided her are available at the Global Biodiversity Information Facility (GBIF, https://doi.org/10.15468/p8gcpe). Data from Michel et al. (2019) are also available at the GBIF (https://doi.org/10.15468/wgfw0h). Stable isotope data from the Belgica 121 and environmental data will be available soon.

Three R scripts are provided:

LinearModels.R is the script assessing the relationships of stable isotope values with trophic groups (factor) and environmental variables (covariables), as well as their first order interactions.

Interactions.R is the script focusing on the relationships of stable isotope values with the first order interactions between continuous environmental covariables.

IsotopicNiches.R is the script investigating the isotopic niches of sea stars depending on trophic groups and environmental variables, through the generation of standard ellipse and the use of original ellipse metrics.
