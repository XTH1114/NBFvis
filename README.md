# NBFvis
A R package for Neighborhood-Based Featurization visualizaion by Tinghui Xu and Kris Sankaran.

## Introduction

Dimensionality reduction of spatial omic data can reveal shared, spatially structured patterns of expression across a collection of genomic features. We study strategies for discovering and interactively visualizing low-dimensional structure in spatial omic data based on the construction of neighborhood features. We design quantile and network-based spatial features that result in spatially consistent embeddings. A simulation compares embeddings made with and without neighborhood-based featurization, and a re-analysis of [Keren et al., 2019] illustrates the overall workflow. We provide an R package, NBFvis, to support computation and interactive visualization for the proposed dimensionality reduction approach.

## Installation

> devtools::install_github("XTH1114/NBFvis")

## Example

Here is the example that reproduces the analysis in our paper.

We first load the packages and dataset we need. The dataset patient4 is a 6643 × 59 data frame of all cells in the tissue section of Patient 4 in the TNBC data [Keren et al., 2019]. We have added two columns named x center and y center, which are the coordinates of the calculated cell centers from the spatial raster data.

> library(NBFvis)
library(dplyr)
data(patient4)

We select 41 variables from dsDNA to HLA Class 1, most of which are proteins and cell type markers. The quantile matrix function generates the quantile matrix from each cell’s neighborhood.

> Quantiles_patient4 <- quantiles_matrix(
  data = patient4 %>% select(dsDNA:HLA_Class_1),
  coordinate = patient4 %>% select(x_center,y_center),
  index = patient4$index,
  NN = 40,
  distance = 60,
  min_percentile = 0.1,
  max_percentile = 0.9,
  quantile_number = 17,
  method = pca_)
  
The function network matrix first builds the network inside the neighborhood and then calculates the corresponding network statistics using the argument given by fun. In this example, we use the function centralities, also exported by our package.

> centrality_patient4 <- network_matrix(
  coordinate = patient4 %>% select(ends_with("_center")),
  index = patient4$index,
  radius = 60,
  NN = 40,
  edge = 30,
  fun = centralities,
  length_output = 29,
  name_output = NULL)

The scales of these two matrix are not the same, which means rescaling is needed. Here we remove Column, index and n neighborhood in the quantile matrix so that all the columns left are quantile and network variables. Normalization and centering are applied to the centralities matrix so that they have a similar scale to the quantile matrix. We then combine the quantile matrix and the rescaled network matrix to construct an extended featurization matrix, which we called the neighborhood matrix earlier.

> neighborhood_info_patient4 <- cbind(
  Quantiles_patient4 %>% select(-index, -n_neighbor),
  scale(centrality_patient4 %>% select(-index)))

The final step is to input the neighborhood matrix, the cell dataset patient4, and the names of the variable of interest in the function NBFvis. This returns an interactive Shiny app that was the source of figures in Section 5.

> NBF_vis(
  matrix = neighborhood_info_patient4,
  origin_data = patient4,
  var_names = colnames(patient4)[17:57])
