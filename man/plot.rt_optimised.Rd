% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oop_conv.R
\name{plot.rt_optimised}
\alias{plot.rt_optimised}
\title{Essentially just does squire::plot_pmcmc_sample. Not intended for use otherwise.}
\usage{
\method{plot}{rt_optimised}(
  x,
  q = c(0.025, 0.975),
  replicates = TRUE,
  summarise = FALSE,
  ci = TRUE,
  particle_fit = FALSE,
  ...
)
}
\arguments{
\item{x}{rt_optimised output object}

\item{q}{Quantiles for upper and lower of interval ribbon}

\item{replicates}{Plot replicates}

\item{summarise}{Logical, add summary line}

\item{ci}{logical add confidence interval ribbon}

\item{particle_fit}{For compatibility, if TRUE deaths are cumulative, else deaths are daily}

\item{...}{placeholder for compatibility does nothing.}
}
\description{
Essentially just does squire::plot_pmcmc_sample. Not intended for use otherwise.
}
