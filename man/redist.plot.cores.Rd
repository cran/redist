% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_cores.R
\name{redist.plot.cores}
\alias{redist.plot.cores}
\title{Plot Cores}
\usage{
redist.plot.cores(shp, plan = NULL, core = NULL, lwd = 2)
}
\arguments{
\item{shp}{A SpatialPolygonsDataFrame or sf object. Required.}

\item{plan}{A numeric vector with one entry for each precinct in shp.
Used to color the districts. Required.}

\item{core}{Required. integer vector produced by \code{redist.identify.cores()}.}

\item{lwd}{Line width. Defaults to 2.}
}
\value{
ggplot
}
\description{
Plot Cores
}
\concept{plot}
