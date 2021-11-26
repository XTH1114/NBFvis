#' Build Networks from Coordinates
#' 
#' This function is used to derive a network object from a coordinate matrix.
#'
#' @param net a data frame containing columns named as \strong{x_center} and \strong{y_center}, which are the x and y coordinates of cells, respectively.
#' @param center_index the index or name of the center cell.
#' @param edge the maximal length for an edge to exist.
#' @return a network object.
#' @importFrom igraph components delete_vertices graph_from_data_frame V
#' @importFrom magrittr %>%
#' @export

build_network <- function(net,center_index,edge){
  n <- nrow(net)
  edgelist <- grid2edgelist(net,edge)
  network1 <- graph_from_data_frame(d = edgelist[,1:2], directed = F)
  sub_gs <- components(network1)$membership
  main_net <- sub_gs[which(names(sub_gs) == as.character(center_index))]
  rm_nodes <- names(which(sub_gs != main_net))
  network1 <- delete_vertices(network1, rm_nodes)
  return(network1)
}
