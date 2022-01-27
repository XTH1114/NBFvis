#' Quantile Matrix Function
#' 
#' This function is used for get the neighborhood quantile matrix for the raw spatial omics dataset. \cr
#' 
#' The spatial omics dataset could be divided into two parts, the expression matrix and the coordinates. \cr
#' Firstly, the dimension reduction method like Principle Component Analysis will be applied to the expression matrix. \cr
#' Then, the neighborhood for each cell will be extracted according to the distance between them. \cr
#' Next, in each neighborhood, tens of quantiles of each principal component are calculated, which are treated as \cr
#' neighborhood information vector for each cell. \cr
#' The final output is a neighborhood quantile matrix, each row of which is the neighborhood quantile vector of \cr
#' principal components for each cell.
#' 
#' 
#' @param data a data frame containing the expression information for each cell, like proteins, genes. \cr
#' This data frame doesn't contains any other columns except expressions, like indices of cells or coordinates of cells.
#' @param coordinate the coordinates of the corresponding cells in the \strong{data} argument. \cr
#' It includes two columns with specific names, the first of which is x coordinates of cells,\cr
#'  and the second of which is the y coordinates of cells.
#' @param index the indices of the corresponding cells in the \strong{data} argument.
#' @param distance the maximal radius of the neighborbood from the center cell. \cr
#' When we try to find a neighborhood for a particular cell, this cell is called the center cell.\cr
#' The cells that are more than the given raidus from the center cell are not be included in the center cell's neighborhood. \cr
#' At least one of \strong{NN} and \strong{dustance} arguments should be specified.
#' @param NN the number of the nearest neighbors. The default value is NULL. \cr
#' It the value is not NULL, the nearest neighborhood method will be applied when we determine the neighborhood for each cell \cr
#' Only the given number of the closest cells will be included in the center cell's neighborhood. \cr
#' At least one of \strong{NN} and \strong{dustance} arguments should be specified.
#' @param min_percentile the minimal percentile. The default value is 0.1.
#' @param max_percentile the maximal percentile. The default value is 0.9.
#' @param quantile_number the number of quantiles for each variable. \cr
#' It is used in the function \strong{seq(min_percentile, max_percentile, length.out = quantile_numbers)} to get the given number of quantiles.
#' @param method the dimension reduction function for the \strong{data} argument. \cr
#' The default method is Principal Component Analysis using \strong{pca_} function in this package.
#' @param ... other parameters passed to the \strong{method} argument.
#' @return a neighborhood quantile matrix. Each row represents quantiles of reduced features in the neighborhood of one cell. \cr
#' The first is the number of cells in the corresponding neighborhood, which is used to check whether there are some abnormal neighborhoods. 
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate filter bind_rows select everything
#' @importFrom rlang set_names
#' @export

quantiles_matrix <- function(data,
                             coordinate,
                             index,
                             distance = NULL,
                             NN = NULL,
                             min_percentile = 0.1,
                             max_percentile = 0.9,
                             quantile_number = 21,
                             method = pca_, ...){
  data <- data[,(sapply(data,var) != 0)]
  tidied_data <- method(data, ...)
  max_column <- ncol(tidied_data)
  colnames(tidied_data) <- paste0("Var",1:max_column)
  n <- nrow(tidied_data)

  final_data <- cbind(x_center = coordinate[,1],
                      y_center = coordinate[,2],
                      tidied_data) %>%
    as.data.frame()

  neighbour <- vector(mode = "list", length = n)
  cat("Generating a list of neighborhoods: \n| ")
  for(i in 1:n){
    if(i %% ceiling(n*0.1) == 0){cat("- ")}
    # xc = x center of the target cell
    xc <- coordinate[i,1]
    yc <- coordinate[i,2]
    temp_1<- final_data%>%
      mutate( dist = sqrt((x_center-xc)^2 + (y_center-yc)^2))
    # head(temp_1)
    if(is.null(distance) == F){
      temp_1 <- temp_1 %>% filter(dist <= distance)
    }
    if(is.null(NN) == F){
      temp_1 <- temp_1[order(temp_1$dist,decreasing = T),]%>%
        .[1:min(NN,nrow(temp_1)),]
    }
    ix <- paste0("Var",1:max_column)
    temp_list <- temp_1[ix]
    neighbour[[i]] <- temp_list
  }
  cat("- |\n")

  cat("Calculating quantiles of each neighborhood: \n| ")
  quantile_list <- vector(mode = "list", length = n)
  for(i in 1:n){
    temp_quantile <- neighbour[[i]] %>%
      apply(2, quantile,probs = seq(min_percentile,max_percentile,length.out = quantile_number)) %>% 
      c()
    temp_quantile <- c(nrow(neighbour[[i]]), temp_quantile) %>%
      set_names(
        c("n_neighbor",
          paste0(
            rep(paste0("Var",1:max_column,"_"),
                each = quantile_number),
            rep(seq(min_percentile,max_percentile,
                    length.out = quantile_number),
                max_column)
          )
        )
      )
    quantile_list[[i]] <- temp_quantile
    if(i %% ceiling(n*0.1) == 0){cat("- ")}
  }
  cat("- |")

  Quantiles <- bind_rows(quantile_list)
  Quantiles$index <- index
  Quantiles <- Quantiles %>% select(index,n_neighbor,everything())

  return(Quantiles)

}

# data <- patient4 %>% select(Na:Au)
# coordinate <- patient4 %>% select(ends_with("_center"))
# quantile_number <- 21
# Quantiles <- quantiles_matrix(data, coordinate,
#                               NN = 40, distance = 100,
#                               index = 39:(38+nrow(data)))



