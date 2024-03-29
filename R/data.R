#' Nimue fit to reported deaths for Afghanistan
#'
#' Included for use in examples. Fitted to data up to 2021-10-06.
#'
#' @format Nimue Simulation Object:
#' \describe{
#'   \item{parameters}{list of fixed parameters}
#'   \item{model}{the odin model}
#'   \item{odin_parameters}{list of fixed parameters in the format needed for the odin model}
#'   \item{replicate_parameters}{data-frame of a sample of the parameters estimated by the MCMC}
#'   \item{pmcmc_results}{full matrix of parameter estimates generated by the pmcmc}
#'   \item{interventions}{various other parameters used in the model, i.e. to generate Rt over-time from the replicates}
#' }
"afg_fit"
#' Nimue fit to excess mortality for Bhutan
#'
#' Included for use in examples. Fitted to data up to 2021-11-16.
#'
#' @format Nimue Simulation Object:
#' \describe{
#'   \item{parameters}{list of fixed parameters}
#'   \item{model}{the odin model}
#'   \item{odin_parameters}{list of fixed parameters in the format needed for the odin model}
#'   \item{replicate_parameters}{data-frame of a sample of the parameters estimated by the MCMC}
#'   \item{pmcmc_results}{full matrix of parameter estimates generated by the pmcmc}
#'   \item{interventions}{various other parameters used in the model, i.e. to generate Rt over-time from the replicates}
#' }
"bhu_fit"
