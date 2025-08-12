# Spatio-temporal-analysis-of-four-decades-of-ant-biodiversity-in-Limburg

This repository contains the statistical method, that is the spatio-temporal model and data cleaning process used for the master thesis "Spatio-temporal analysis of four decades of ant biodiversity in Limburg". Due to data confidentiality the data provided in this github repository differs from the one used in the master thesis.

## Different files
This repository contains three different R files, which should be applied to the data in the following order:
1) data_cleaning.R
2) data_exploration.R
3) single_species_analysis.R

An explanation on how to use each R file is given below. However first some general information and recommendations are given.

## General set up
Each R file consists of a set up section, which consists of general parameters that can be varied according to the users preference or expert/historical knowledge. In this set up section, the required packages are loaded first (note that these need to be installed before hand. For user using a Unix operating systems, note that dependencies on other R packages must be specified when installing an R package). Secondly, data and the user specified parameters are loaded (Here you need to specify certain parameters, otherwise the code would not work). Third a file specific section is given which contains the data cleaning (to satisfy the required inclusion criteria, that are defined by the user), a global data exploration, single species analysis respectively, for the denoted R files above. 

## R-file specific information
### data_cleaning.R usage

### data_exploration.R usage

### single_species_analysis.R
