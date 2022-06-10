library(NBFvis)
library(dplyr)
data(patient4)
Quantiles_patient4 <- quantiles_matrix(
  data = patient4 %>% select(dsDNA:HLA_Class_1),
  coordinate = patient4 %>% select(x_center,y_center),
  index = patient4$index,
  NN = 40,
  distance = 60,
  min_percentile = 0.1,
  max_percentile = 0.9,
  quantile_number = 17,
  method = pca_)

centrality_patient4 <- network_matrix(
  coordinate = patient4 %>% select(ends_with("_center")),
  index = patient4$index,
  radius = 60,
  NN = 40,
  edge = 30,
  fun = centralities,
  length_output = 29,
  name_output = NULL)

neighborhood_info_patient4 <- cbind(
  Quantiles_patient4 %>% select(-index, -n_neighbor),
  scale(centrality_patient4 %>% select(-index)))

NBF_vis(
  matrix = neighborhood_info_patient4,
  origin_data = patient4,
  var_names = colnames(patient4)[17:57])
