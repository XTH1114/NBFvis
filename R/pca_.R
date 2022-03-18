#' Principal Component Analysis function used in the function quantile_matrix 
#'
#' This function is used for Principal Component Analysis(PCA) in the \cr
#' function \strong{quantile_matrix} in this package. It uses \strong{prcomp} function for PCA.
#'
#' @param data the input data matrix, containing all the expression information like proteins or genes. \cr
#' Other things like indices and coordinates should not be included.
#' @param pca_threshold the minimal percentage of the explained variance in Principal Component Analysis. \cr
#' The default value is 0.9.
#' @param ... other parameters could be passed to the \strong{prcomp} function.
#' @export

pca_ <- function(data, pca_threshold = 0.9, ...){
  data_pca <- prcomp(data, ...)
  max_column <- min(which(summary(data_pca)$importance[3,]>pca_threshold))
  tidied_data <- data_pca$x[,1:max_column]
  return(tidied_data)
}
