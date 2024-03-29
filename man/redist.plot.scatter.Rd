% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_plans.R
\name{redist.plot.scatter}
\alias{redist.plot.scatter}
\title{Scatter plot of plan summary statistics}
\usage{
redist.plot.scatter(plans, x, y, ..., bigger = TRUE)
}
\arguments{
\item{plans}{the \code{redist_plans} object.}

\item{x}{\code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the
quantity to plot on the horizontal axis.}

\item{y}{\code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the
quantity to plot on the vertical axis.}

\item{...}{passed on to \code{\link[ggplot2:geom_point]{geom_point}}.}

\item{bigger}{if TRUE, make the point corresponding to the reference plan larger.}
}
\value{
A ggplot
}
\description{
Makes a scatterplot of two quantities of interest across districts or plans.
}
\examples{
library(dplyr)
data(iowa)

iowa <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05, total_pop = pop)
plans <- redist_smc(iowa, nsims = 100, silent = TRUE)
plans \%>\%
    mutate(comp = distr_compactness(iowa)) \%>\%
    group_by(draw) \%>\%
    summarize(pop_dev = max(abs(total_pop/mean(total_pop) - 1)),
        comp = comp[1]) \%>\%
    redist.plot.scatter(pop_dev, comp)

}
\concept{plot}
