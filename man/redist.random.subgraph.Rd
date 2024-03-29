% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random_subgraph.R
\name{redist.random.subgraph}
\alias{redist.random.subgraph}
\title{Return a random subgraph of a shape}
\usage{
redist.random.subgraph(shp, n, adj = NULL)
}
\arguments{
\item{shp}{sf object or SpatialPolygonsDataFrame}

\item{n}{number of edges to sample. n must be a positive integer.}

\item{adj}{Optional. zero indexed adjacency list.}
}
\value{
sf dataframe with n rows
}
\description{
\code{random.subgraph} returns a random subset of the shp provided
}
\details{
Snowball sampling with backtracking
}
