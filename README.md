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

In order to smoothly run the code, the specified hardware requirements below are advised. Especially the RAM storage becomes important, once the model validation starts. 
Hardware requirements: 
1) minimum of 32 GB RAM
2) minimum of CPU with 4 cores 

## R-file specific information
### data_cleaning.R usage
This R-file contains the Rcode to properly set up the map, identify the different habitats of all observed points, add presence/absence variable, seasonality effect, etc. This code must be run first in order to continou with the data exploration and data analysis. Note that due to confidentiallity this R-file is highly restrictive to the used dataset (which is not the dataset that is available in this repository). The dataset in this repository consists of public available data. This dataset was already modified and thus this R-file can be neglegted. This R-file also implements the inclusion criteria to the data, such that the data provided to us is ready for data exploration and analysis.

### data_exploration.R usage
This R-file contains the Rcode to give a general idea of weaknesses and strengths of the data, by exploring various variables, such as habitat type distribution, endangered species distibution, seasonality distribution, etc. .

### single_species_analysis.R
This R-file contains the Rcode to perfomr the statistical analysis described in the Master thesis. Single species occupancy models were constructed by the use of the INLA package. Various models were fit to the modified data. This models all consists of the following:
1) Data preparation (modifying the presence/absence variable such that this is species specific).
2) Mesh construction
3) SPDE approach to represent GRF with Matèrn covariance structure as a GMRF with Matèrn covariance structure.
4) Index set to link each mesh node with each time point. It thus provide structure to the spatio-temporal random effect.
5) Projection matrix to project the continouos GRF with Matèrn covariacne to the mesh nodes
6) Prediction data that constructs locations and times where predictions must be made
7) Stacking all previous steps, such that INLA can eficiently select all the created components
8) Model formulation
9) Model fitting using INLA
