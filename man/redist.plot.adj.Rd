% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_map.R
\name{redist.plot.adj}
\alias{redist.plot.adj}
\title{Creates a Graph Overlay}
\usage{
redist.plot.adj(
  shp,
  adj = NULL,
  plan = NULL,
  centroids = TRUE,
  drop = FALSE,
  plot_shp = TRUE,
  zoom_to = NULL,
  title = ""
)
}
\arguments{
\item{shp}{A SpatialPolygonsDataFrame or sf object. Required.}

\item{adj}{A zero-indexed adjacency list. Created with redist.adjacency
if not supplied. Default is NULL.}

\item{plan}{A numeric vector with one entry for each precinct in shp.
Used to remove edges that cross boundaries. Default is \code{NULL}.  Optional.}

\item{centroids}{A logical indicating if centroids should be plotted. Default is \code{TRUE}.}

\item{drop}{A logical indicating if edges that cross districts should be dropped. Default is \code{FALSE}.}

\item{plot_shp}{A logical indicating if the shp should be plotted under the
graph. Default is \code{TRUE}.}

\item{zoom_to}{\code{\link[dplyr:dplyr_data_masking]{<data-masking>}} An
indexing vector of units to zoom the map to.}

\item{title}{A string title of plot. Defaults to empty string. Optional.}
}
\value{
ggplot map
}
\description{
Creates a Graph Overlay
}
\examples{
data(iowa)
redist.plot.adj(shp = iowa, plan = iowa$cd_2010)

}
\concept{plot}
