% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmcmc_likelihood_excess.R
\name{excess_log_likelihood}
\alias{excess_log_likelihood}
\title{Likelihood function for excess mortality data
Not for individual use, only exported so it is easy to call in pmcmc functions.}
\usage{
excess_log_likelihood(
  pars,
  data,
  squire_model,
  model_params,
  pars_obs,
  n_particles,
  forecast_days = 0,
  return = "ll",
  Rt_args,
  interventions,
  ...
)
}
\arguments{
\item{pars}{Parameter inputs}

\item{data}{data to fit to}

\item{squire_model}{the model object}

\item{model_params}{the fixed parameters}

\item{pars_obs}{extra non fitted parameters}

\item{n_particles}{unused}

\item{forecast_days}{unused}

\item{return}{Specifies we are return likelihood for model outputs}

\item{Rt_args}{parameters for Rt}

\item{interventions}{parameters around vaccination}

\item{...}{unused (for compatibility)}
}
\description{
Likelihood function for excess mortality data
Not for individual use, only exported so it is easy to call in pmcmc functions.
}
