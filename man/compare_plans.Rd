% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify_compare.R
\name{compare_plans}
\alias{compare_plans}
\title{Make a comparison between two sets of plans}
\usage{
compare_plans(
  plans,
  set1,
  set2,
  shp = NULL,
  plot = "fill",
  thresh = 0.1,
  labs = c("Set 1", "Set 2"),
  ncores = 1
)
}
\arguments{
\item{plans}{a \link{redist_plans} object}

\item{set1}{\code{\link[dplyr:dplyr_data_masking]{<data-masking>}} indexing vectors
for the plan draws to compare.  Alternatively, a second \link{redist_plans}
object to compare to.}

\item{set2}{\code{\link[dplyr:dplyr_data_masking]{<data-masking>}} indexing vectors
for the plan draws to compare. Must be mutually exclusive with \code{set1}.}

\item{shp}{a shapefile for plotting.}

\item{plot}{If \code{plot="line"}, display a plot for each set showing the set of
boundaries which most distinguish it from the other set (the squared
differences in the eigenvector values across the boundary).  If
\code{plot="fill"}, plot the eigenvector for each set as a choropleth. If \code{plot = 'adj'},
plot the shows the adjacency graph edges which most distinguish it from the other set.
The adj option is a different graphical option of the same information as the line
option. See below for more information.  Set to \code{FALSE} to disable plotting
(or leave out \code{shp}).}

\item{thresh}{the value to threshold the eigenvector at in determining the
relevant set of precincts for comparison.}

\item{labs}{the names of the panels in the plot.}

\item{ncores}{the number of parallel cores to use.}
}
\value{
If possible, makes a comparison plot according to \code{plot}. Otherwise
returns the following list:
\item{eigen1}{A numeric vector containing the first eigenvector of
\code{p1 - p2}, where \code{p1} and \code{p2} are the co-occurrence matrices
for \code{set1} and \code{set2}, respectively.}
\item{eigen2}{A numeric vector containing the first eigenvector of
\code{p2 - p1}, where \code{p1} and \code{p2} are the co-occurrence matrices
for \code{set1} and \code{set2}, respectively.}
\item{group_1a, group_1b}{Lists of precincts. Compared to \code{set2}, in the
\code{set1} plans these precincts were much more likely to be in separate
districts. Computed by thresholding \code{eigen1} at \code{thresh}.}
\item{group_2a, group_2b}{Lists of precincts. Compared to \code{set1}, in the
\code{set2} plans these precincts were much more likely to be in separate
districts. Computed by thresholding \code{eigen2} at \code{thresh}.}
\item{cooccur_sep_1}{The difference in the average co-occurrence of precincts
in \code{group_1a} and \code{group_1b} between \code{set2} and \code{set1}.
Higher indicates better separation.}
\item{cooccur_sep_2}{The difference in the average co-occurrence of precincts
in \code{group_2a} and \code{group_2b} between \code{set1} and \code{set2}.
Higher indicates better separation.}
}
\description{
This function provides one way to identify the structural differences between
two sets of redistricting plans. It operates by computing the precinct
co-occurrence matrix (a symmetric matrix where the i,j-th entry is the
fraction of plans where precinct i and j are in the same district) for each
set, and then computing the first eigenvalue of the difference in these two
matrices (in each direction). These eigenvalues identify the important parts
of the map.
}
\details{
The co-occurrence matrices are regularized with a \eqn{Beta(1/ndists, 1-1/ndists)}
prior, which is useful for when either \code{set1} or \code{set2} is small.
}
\examples{
data(iowa)
iowa_map <- redist_map(iowa, ndists = 4, pop_tol = 0.05)
plans1 <- redist_smc(iowa_map, 100, silent = TRUE)
plans2 <- redist_mergesplit(iowa_map, 200, warmup = 100, silent = TRUE)
compare_plans(plans1, plans2, shp = iowa_map)
compare_plans(plans2, as.integer(draw) <= 20,
    as.integer(draw) > 20, shp = iowa_map, plot = "line")

}
\concept{analyze}
