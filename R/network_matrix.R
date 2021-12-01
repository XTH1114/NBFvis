#' Calculate Network Statistics Matrix
#'
#' This function is aimed at calculating network statistics of the neighborhood \cr
#' network of each cell. The output is a network statistics matrix.
#'
#' @param coordinate a data frame containing cell indices and x and y coordinates which are labeled as \strong{index}, \strong{x_center} and \strong{y_center}, respectively.
#' @param radius the maximal radius length of a neighborhood. The cells whose distance to the center cell is larger than this length are not \cr
#' included in the neighborhood.
#' @param NN the number of nearest neighbors. If not NULL, only the nearest NN cells will be included in neighborhood.
#' @param edge the maximal length of the network edge. 
#' @param fun the function which calculates the statistics of a network. \cr
#' Two arguments are supposed to be in this function. The first one is a network object \cr
#' and the second one is an index indicating the statistics of which vertex (cell) should be calculated. \cr
#' The default function is \strong{centralities}.
#' @param length_output the number of statistics calculated in the argument \strong{fun}. \cr
#' The default value is 30, which is the number of statistics in the \strong{centralities} function of this package.
#' @param name_output the names of statistics calculated in the argument \strong{fun}.
#' @return a network statistics matrix, each rows of it representing network statistics of the neighborhood network of a single cell.
#' @importFrom dplyr mutate select everything filter bind_rows
#' @importFrom igraph V
#' @importFrom purrr set_names
#' @import expm
#' @export
#'



network_matrix <- function(coordinate,
                           index,
                           radius = NULL,
                           NN = NULL,
                           edge,
                           fun = centralities,
                           length_output = 30,
                           name_output = NULL){

  n <- nrow(coordinate)
  centrality_table <- vector(mode = "list", length = n)
  cat("Calculating network statistics with respect to the center cell: \n| ")
  perc <- 0
  for (i in 1:n){
    center_cell <- coordinate[i,]
    xc <- center_cell$x_center
    yc <- center_cell$y_center
    net1 <- coordinate%>%
      mutate(distance = sqrt((x_center-xc)^2+(y_center-yc)^2))
    net1 <- net1 %>% mutate(id = index) %>%
      select(id, x_center,y_center,distance,everything())
    if(is.null(radius) == F){
      net1 <- net1 %>% filter(distance<=radius)
    }
    if(is.null(NN) == F){
      net1 <- net1[which((rank(net1$distance))<=(NN+1)),]}

    network1 <- build_network(net1,index[i],edge)

    if(i %% floor(length(index)*0.01) == 0){
      if(perc< 100){cat("- ")}
      if((perc+1) %% 10 == 0 & perc != 0){cat(paste0((perc+1) %/% 10,"0%"))}
      if((perc+1) %% 20 == 0 & (perc+1) != 100 & (perc+1) != 0){cat(" |\n| ")}
      if((perc+1) == 100){cat("|")}
      perc <- perc+1}

    if(length(V(network1)) == 1){
      centrality_table[[i]] <- c(rep(NA,times = length_output))
    } else{
      output <- fun(network1,center_index = as.character(index[i]))
      if(is.null(name_output) == F){
        centrality_table[[i]] <- output %>% set_names(nm = name_output)
      }else{
        centrality_table[[i]] <- output %>% set_names(nm = paste0("statistics_", 1:length_output))
      }
    }
  }
  return(cbind(index = index, bind_rows(centrality_table)))

}

# load("patient4.Rda")
# coordinate <- patient4 %>% select(ends_with("_center"))
# index <- 39:(nrow(coordinate)+38)
# radius = 100
# NN = 40
# edge = 40
# name_output <-  NULL
# length_output <- 30
#
# c <- network_matrix(coordinate,index,radius = 100,NN = 40,edge,
#                     fun = centralities,length_output = 30,name_output = NULL)
#
