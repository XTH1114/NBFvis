% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grid2edgelist.R
\name{grid2edgelist}
\alias{grid2edgelist}
\title{Get the Edgelist from Coordinates}
\usage{
grid2edgelist(grid, edge)
}
\arguments{
\item{grid}{the coordinates of cells. It includes at least three columns, \strong{x_center}, \strong{y_center}, and \strong{id}. \cr
\strong{x_center} is the x coordinates of the cells, and \strong{y_center} is the y_coordinates of the cells. \cr
\strong{id} is the indices of thec cells.}

\item{edge}{the maximal length of an edge in the network.}
}
\value{
an edgelist of a network generated from coordinates of cells using the given maximal edge length.
}
\description{
This function is used to get the edgelist of a network from coordinates of cells. \cr
A network is built according to the coordinates of cells. \cr
Edges will be generated between the nodes whose distance is smaller than or equal to the given length. \cr
This function is used in the function \strong{build_network}
}
\details{
The output is an edgelist of this network.
}
