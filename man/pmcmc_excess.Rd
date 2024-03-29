% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmcmc_excess.R
\name{pmcmc_excess}
\alias{pmcmc_excess}
\title{Runs an adaptive metropolis hasting MCMC on the specified model.
Adapted from squire::pmcmc, but handles the weekly death data from excess-mortality.}
\usage{
pmcmc_excess(
  data,
  n_mcmc,
  log_likelihood = NULL,
  log_prior = NULL,
  use_drjacoby = FALSE,
  drjacoby_list = NULL,
  n_particles = 100,
  steps_per_day = 1,
  output_proposals = FALSE,
  n_chains = 1,
  squire_model = squire::explicit_model(),
  pars_obs = list(phi_cases = 1, k_cases = 2, phi_death = 1, k_death = 2, exp_noise =
    1e+07, likelihood = function(model_deaths, data_deaths) {
    
    squire:::ll_nbinom(data_deaths, model_deaths, pars_obs$phi_death, pars_obs$k_death,
    pars_obs$exp_noise)
 }),
  pars_init = list(start_date = as.Date("2020-02-07"), R0 = 2.5, Meff = 2, Meff_pl = 3,
    R0_pl_shift = 0),
  pars_min = list(start_date = as.Date("2020-02-01"), R0 = 0, Meff = 1, Meff_pl = 2,
    R0_pl_shift = -2),
  pars_max = list(start_date = as.Date("2020-02-20"), R0 = 5, Meff = 3, Meff_pl = 4,
    R0_pl_shift = 5),
  pars_discrete = list(start_date = TRUE, R0 = FALSE, Meff = FALSE, Meff_pl = FALSE,
    R0_pl_shift = FALSE),
  proposal_kernel = NULL,
  scaling_factor = 1,
  reporting_fraction = 1,
  country = NULL,
  population = NULL,
  contact_matrix_set = NULL,
  baseline_contact_matrix = NULL,
  date_contact_matrix_set_change = NULL,
  date_Rt_change = NULL,
  hosp_bed_capacity = NULL,
  baseline_hosp_bed_capacity = NULL,
  date_hosp_bed_capacity_change = NULL,
  ICU_bed_capacity = NULL,
  baseline_ICU_bed_capacity = NULL,
  date_ICU_bed_capacity_change = NULL,
  date_vaccine_change = NULL,
  baseline_max_vaccine = NULL,
  max_vaccine = NULL,
  date_vaccine_efficacy_infection_change = NULL,
  baseline_vaccine_efficacy_infection = NULL,
  vaccine_efficacy_infection = NULL,
  date_vaccine_efficacy_disease_change = NULL,
  baseline_vaccine_efficacy_disease = NULL,
  vaccine_efficacy_disease = NULL,
  Rt_args = NULL,
  burnin = 0,
  replicates = 100,
  forecast = 0,
  required_acceptance_ratio = 0.23,
  start_adaptation = round(n_mcmc/2),
  gibbs_sampling = FALSE,
  gibbs_days = NULL,
  dur_R = 365,
  dose_ratio = NULL,
  date_vaccine_efficacy = NULL,
  vaccine_efficacies = NULL,
  ...
)
}
\arguments{
\item{data}{Data to fit to.  This must be constructed with
\code{particle_filter_data}}

\item{n_mcmc}{number of mcmc mcmc iterations to perform}

\item{log_likelihood}{function to calculate log likelihood, must take named parameter vector as input,
allow passing of implicit arguments corresponding to the main function arguments.
Returns a named list, with entries:
  - $log_likelihood, a single numeric
  - $sample_state, a numeric vector corresponding to the state of a single particle, chosen at random,
  at the final time point for which we have data.
  If NULL, calculated using the function calc_loglikelihood.}

\item{log_prior}{function to calculate log prior, must take named parameter vector as input, returns a single numeric.
If NULL, uses uninformative priors which do not affect the posterior}

\item{use_drjacoby}{Should we use drjacoby to fit the model}

\item{drjacoby_list}{Parameters to use in drjacoby}

\item{n_particles}{Number of particles (considered for both the PMCMC fit and sampling from posterior)}

\item{steps_per_day}{Number of steps per day}

\item{output_proposals}{Logical indicating whether proposed parameter jumps should be output along with results}

\item{n_chains}{number of MCMC chains to run}

\item{squire_model}{A squire model to use}

\item{pars_obs}{list of parameters to use in comparison
with \code{compare}. Must be a list containing, e.g.
list(phi_cases = 0.1,
     k_cases = 2,
     phi_death = 1,
     k_death = 2,
     exp_noise = 1e6)}

\item{pars_init}{named list of initial inputs for parameters being sampled}

\item{pars_min}{named list of lower reflecting boundaries for parameter proposals}

\item{pars_max}{named list of upper reflecting boundaries for parameter proposals}

\item{pars_discrete}{named list of logicals, indicating if proposed jump should be discrete}

\item{proposal_kernel}{named matrix of proposal covariance for parameters}

\item{scaling_factor}{numeric for starting scaling factor for covariance matrix. Default = 1}

\item{reporting_fraction}{Reporting fraction. Numeric for what proportion of
the total deaths the reported deaths represent. E.g. 0.5 results in
the model calibrating to twice the deaths provided by \code{data$deaths}}

\item{country}{Character for country beign simulated. WIll be used to
generate \code{population} and \code{contact_matrix_set} if
unprovided. Either \code{country} or \code{population} and
\code{contact_matrix_set} must be provided.}

\item{population}{Population vector (for each age group). Default = NULL,
which will cause population to be sourced from \code{country}}

\item{contact_matrix_set}{List of contact matrices to be used from the dates
provided in \code{date_contact_matrix_set_change}.Must be same length as
\code{date_contact_matrix_set_change}}

\item{baseline_contact_matrix}{The starting contact matrix prior to any changes
due to interventions or otherwise. Default = NULL, which will use the contact
matrix associated with the coutnry provided.}

\item{date_contact_matrix_set_change}{Calendar dates at which the contact matrices
set in \code{model_params} change. Defaut = NULL, i.e. no change}

\item{date_Rt_change}{The dates corresponding to each Rt change}

\item{hosp_bed_capacity}{Number of hospital beds at each date specified in
\code{date_hosp_bed_capacity_change}. Must be same length as
\code{date_hosp_bed_capacity_change}.}

\item{baseline_hosp_bed_capacity}{The starting number of hospital beds before
the epidemic started. Default = NULL, which will use the hospital beds data
for the country provided. If no country is provided then this is 5/1000 of
the population}

\item{date_hosp_bed_capacity_change}{Calendar dates at which hospital bed
capacity changes set in \code{model_params} change.
Defaut = NULL, i.e. no change}

\item{ICU_bed_capacity}{Number of ICU beds at each date specified in
\code{date_ICU_bed_capacity_change}. Must be same length as
\code{date_ICU_bed_capacity_change}.}

\item{baseline_ICU_bed_capacity}{The starting number of ICU beds before
the epidemic started. Default = NULL, which will use the hospital beds data
for the country provided. If no country is provided then this is 3/100 of
hospital beds}

\item{date_ICU_bed_capacity_change}{Calendar dates at which ICU bed
capacity changes set in \code{model_params} change.
Defaut = NULL, i.e. no change}

\item{date_vaccine_change}{Date that vaccine doses per day change.
Default = NULL.}

\item{baseline_max_vaccine}{Baseline vaccine doses per day. Default = NULL}

\item{max_vaccine}{Time varying maximum vaccine doeses per day. Default = NULL.}

\item{date_vaccine_efficacy_infection_change}{Date that vaccine efficacy
against infection changes. Default = NULL.}

\item{baseline_vaccine_efficacy_infection}{Baseline vaccine effacy against infection.
Default = NULL}

\item{vaccine_efficacy_infection}{Time varying vaccine efficacy against infection.
Default = NULL.}

\item{date_vaccine_efficacy_disease_change}{Date that vaccine efficacy
against disease changes. Default = NULL.}

\item{baseline_vaccine_efficacy_disease}{Baseline vaccine efficacy against disease
Default = NULL}

\item{vaccine_efficacy_disease}{Time varying vaccine efficacy against infection.
Default = NULL.}

\item{Rt_args}{List of arguments to be passed to \code{evaluate_Rt_pmcmc} for calculating Rt.
Current arguments are available in \code{Rt_args_list}}

\item{burnin}{number of iterations to discard from the start of MCMC run when sampling from the posterior for trajectories}

\item{replicates}{number of trajectories (replicates) to be returned that are being sampled from the posterior probability results produced by \code{run_mcmc_chain}
to select parameter set. For each parmater set sampled, run particle filter with \code{n_particles} and sample 1 trajectory}

\item{forecast}{Number of days to forecast forward. Default = 0}

\item{required_acceptance_ratio}{Desired MCMC acceptance ratio}

\item{start_adaptation}{Iteration number to begin RM optimisation of scaling factor at}

\item{gibbs_sampling}{Whether or not to use the Gibbs Sampler for start_date}

\item{gibbs_days}{Number of days either side of the start_date parameter to evaluate likelihood at}

\item{dur_R}{Allows the specification of the duration of immunity}

\item{dose_ratio}{The time dependent number of people with two doses}

\item{date_vaccine_efficacy}{Dates at which VE changes}

\item{vaccine_efficacies}{Values for the VE that can vary over time.}

\item{...}{Further aguments for the model parameter function. If using the
\code{\link[squire]{explicit_model}} (default) this will be
\code{parameters_explicit_SEEIR}.}
}
\description{
Runs an adaptive metropolis hasting MCMC on the specified model.
Adapted from squire::pmcmc, but handles the weekly death data from excess-mortality.
}
