#' A Function Calculating Network Centralities Statistics
#'
#'This function is used to calculate the network statistics with respect to one specific vertex \cr
#'in a network, most of which are various centralities.
#'
#' @param g A network object (g stands for "graph").
#' @param center_index The index or name of the vertex (center cell) whose centralities we are interested in.
#' @importFrom igraph V
#' @importFrom igraph degree betweenness closeness eigen_centrality eccentricity subgraph_centrality get.adjacency
#' @importFrom sna loadcent gilschmidt infocent stresscent
#' @importFrom centiserve averagedis barycenter closeness.latora closeness.residual
#' @importFrom centiserve communibet crossclique decay diffusion.degree entropy
#' @importFrom centiserve geokpath laplacian leverage lincent lobby markovcent mnc radiality semilocal topocoefficient
#' @export
#' @references https://www.r-bloggers.com/2018/12/network-centrality-in-r-an-introduction/

centralities <- function(g, center_index){
  num_stats <- 29
  statistics <- vector(mode = "numeric", num_stats)
  if(!(as.character(center_index) %in% V(g)$name)){
    return(rep(NA,num_stats))
  }
  relative_index <- V(g)$name == as.character(center_index)

  #igraph
  statistics[1] <- length(V(g))
  statistics[2] <- degree(g,v = center_index)
  statistics[3] <- betweenness(g,v = center_index)
  statistics[4] <- closeness(g,vids = center_index)
  statistics[5] <- eigen_centrality(g)$vector[relative_index]
  statistics[6] <- 1/eccentricity(g,vids = center_index)
  statistics[7] <- subgraph_centrality(g)[relative_index]

  A <- get.adjacency(g,sparse=F)
 #sna
  statistics[8] <- loadcent(A,nodes = relative_index)
  statistics[9] <- gilschmidt(A,nodes = relative_index)
  statistics[10] <- infocent(A,nodes = relative_index)
  statistics[11] <- stresscent(A,nodes = relative_index)

  # centiserve
  statistics[12] <- 1/averagedis(g,vids = center_index)
  statistics[13] <- barycenter(g,vids = center_index)
  statistics[14] <- closeness.latora(g,vids = center_index)
  statistics[15] <- closeness.residual(g,vids = center_index)
  statistics[16] <- communibet(g,vids = center_index)
  statistics[17] <- crossclique(g,vids = center_index)
  statistics[18] <- decay(g,vids = center_index)
  statistics[19] <- diffusion.degree(g,vids = center_index)
  statistics[20] <- geokpath(g,vids = center_index)
  statistics[21] <- laplacian(g,vids = center_index)
  statistics[22] <- leverage(g,vids = center_index)
  statistics[23] <- lincent(g,vids = center_index)
  statistics[24] <- lobby(g,vids = center_index)
  statistics[25] <- markovcent(g,vids = center_index)
  statistics[26] <- mnc(g,vids = center_index)
  statistics[27] <- radiality(g,vids = center_index)
  statistics[28] <- semilocal(g,vids = center_index)
  statistics[29] <- 1/topocoefficient(g,vids = center_index)
  return(round(statistics,4))
}
