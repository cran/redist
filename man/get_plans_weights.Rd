% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/redist_plans.R
\name{get_plans_weights}
\alias{get_plans_weights}
\alias{weights.redist_plans}
\title{Extract the sampling weights from a redistricting simulation.}
\usage{
get_plans_weights(plans)

\method{weights}{redist_plans}(object, ...)
}
\arguments{
\item{plans, object}{the \code{redist_plans} object}

\item{...}{Ignored.}
}
\value{
A numeric vector of weights, with an additional attribute
\code{resampled} indicating whether the plans have been resampled according
to these weights. If weights have been resampled, this returns the weights
before resampling (i.e., they do not correspond to the resampled plans).

numeric vector
}
\description{
May be \code{NULL} if no weights exist (MCMC or optimization methods).
}
\concept{analyze}
