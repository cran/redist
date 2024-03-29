% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_diagnostics.R
\name{redist.diagplot}
\alias{redist.diagplot}
\title{Diagnostic plotting functionality for MCMC redistricting.}
\usage{
redist.diagplot(sumstat,
plot = c("trace", "autocorr", "densplot", "mean", "gelmanrubin"),
logit = FALSE, savename = NULL)
}
\arguments{
\item{sumstat}{A vector, list, \code{mcmc} or \code{mcmc.list} object
containing a summary statistic of choice.}

\item{plot}{The type of diagnostic plot to generate: one of "trace",
"autocorr", "densplot", "mean", "gelmanrubin". If \code{plot = "gelmanrubin"},
the input \code{sumstat} must be of class \code{mcmc.list} or \code{list}.}

\item{logit}{Flag for whether to apply the logistic transformation for the
summary statistic. The default is \code{FALSE}.}

\item{savename}{Filename to save the plot. Default is \code{NULL}.}
}
\value{
Returns a plot of file type \code{.pdf}.
}
\description{
\code{redist.diagplot} generates several common MCMC diagnostic plots.
}
\details{
This function allows users to generate several standard diagnostic
plots from the MCMC literature, as implemented by Plummer et. al (2006).
Diagnostic plots implemented include trace plots, autocorrelation plots,
density plots, running means, and Gelman-Rubin convergence diagnostics
(Gelman & Rubin 1992).
}
\examples{
\donttest{
data(fl25)
data(fl25_enum)
data(fl25_adj)

## Get an initial partition
init_plan <- fl25_enum$plans[, 5118]
fl25$init_plan <- init_plan

## 25 precinct, three districts - no pop constraint ##
fl_map <- redist_map(fl25, existing_plan = 'init_plan', adj = fl25_adj)
alg_253 <- redist_flip(fl_map, nsims = 10000)

## Get Republican Dissimilarity Index from simulations
rep_dmi_253 <- redistmetrics::seg_dissim(alg_253, fl25, mccain, pop) |>
    redistmetrics::by_plan(ndists = 3)

## Generate diagnostic plots
redist.diagplot(rep_dmi_253, plot = "trace")
redist.diagplot(rep_dmi_253, plot = "autocorr")
redist.diagplot(rep_dmi_253, plot = "densplot")
redist.diagplot(rep_dmi_253, plot = "mean")

## Gelman Rubin needs two chains, so we run a second
alg_253_2 <- redist_flip(fl_map, nsims = 10000)

rep_dmi_253_2 <- redistmetrics::seg_dissim(alg_253_2, fl25, mccain, pop) |>
    redistmetrics::by_plan(ndists = 3)

## Make a list out of the objects:
rep_dmi_253_list <- list(rep_dmi_253, rep_dmi_253_2)

## Generate Gelman Rubin diagnostic plot
redist.diagplot(sumstat = rep_dmi_253_list, plot = "gelmanrubin")

}
}
\references{
Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander
Tarr. (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
Carlo." Working Paper. Available at
\url{http://imai.princeton.edu/research/files/redist.pdf}.

Gelman, Andrew and Donald Rubin. (1992) "Inference from iterative simulations
using multiple sequences (with discussion)." Statistical Science.

Plummer, Martin, Nicky Best, Kate Cowles and Karen Vines. (2006) "CODA:
Convergence Diagnosis and Output Analysis for MCMC." R News.
}
\concept{plot}
