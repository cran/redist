% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enumpart.R
\name{redist.read.enumpart}
\alias{redist.read.enumpart}
\title{Read Results from enumpart}
\usage{
redist.read.enumpart(out_path, skip = 0, n_max = -1L)
}
\arguments{
\item{out_path}{out_path specified in redist.run.enumpart}

\item{skip}{number of lines to skip}

\item{n_max}{max number of lines to read}
}
\value{
district_membership matrix
}
\description{
Read Results from enumpart
}
\examples{
\dontrun{
temp <- tempdir()
cds <- redist.read.enumpart(out_path = paste0(temp, "/enumerated"))
}
}
\references{
Benjamin Fifield, Kosuke Imai, Jun Kawahara, and Christopher T Kenny.
"The Essential Role of Empirical Validation in Legislative Redistricting Simulation."
Forthcoming, Statistics and Public Policy.
}
\concept{enumerate}
