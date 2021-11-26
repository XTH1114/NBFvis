#' Get the Edgelist from Coordinates 
#' 
#' This function is used to get the edgelist of a network from coordinates of cells. \cr
#' A network is built according to the coordinates of cells. \cr
#' Edges will be generated between the nodes whose distance is smaller than or equal to the given length. \cr
#' This function is used in the function \strong{build_network}
#' 
#' The output is an edgelist of this network.
#' @param grid the coordinates of cells. It includes at least three columns, \strong{x_center}, \strong{y_center}, and \strong{id}. \cr
#' \strong{x_center} is the x coordinates of the cells, and \strong{y_center} is the y_coordinates of the cells. \cr
#' \strong{id} is the indices of thec cells.
#' @param edge the maximal length of an edge in the network.
#' @return an edgelist of a network generated from coordinates of cells using the given maximal edge length.
#' @importFrom magrittr %>% 
#' @importFrom igraph V
#' @importFrom purrr set_names
#' @importFrom dplyr bind_rows
#' @export
#' 


grid2edgelist <- function(grid, edge){
  n <- nrow(grid)
  id <- grid$id
  temp_list <- vector(mode = "list", length = n-1)
  for(i in 1:(n-1)){
    node1 <- c(grid$x_center[i],grid$y_center[i])
    dists <- sqrt((grid[(i+1):n,c("x_center")] - node1[1])^2+(grid[(i+1):n,c("y_center")] - node1[2])^2)
    dists_exist_index <- which(dists <= edge)
    dists_left_index <- id[(i+1):n][dists_exist_index]
    dists_left <- dists[dists_exist_index]
    temp_list[[i]] <- data.frame(rep(id[i],length(dists_left)),
                            dists_left_index,
                            dists_left) %>%
      set_names(nm = c("edge1", "edge2", "distance"))
  }
  return(bind_rows(temp_list))
}
