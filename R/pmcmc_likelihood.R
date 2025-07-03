#' Methods to get the log-likelihood of a model
#' @param model_out Model object
#' @noRd
get_model_likelihood <- function(model_out){
  UseMethod("get_model_likelihood")
}
#' Method to calculate log-likelihood for each model type for a nimue model
#' @param model_out Model object
#' @export
get_model_likelihood.nimue_simulation <- function(model_out){
  squire:::calc_loglikelihood
}
#' Method to calculate log-likelihood for each model type for a nimue model using excess deaths
#' @param model_out Model object
#' @export
get_model_likelihood.excess_nimue_simulation <- function(model_out){
  excess_log_likelihood
}
#' Method to calculate log-likelihood for each model type for a nimue model using excess deaths with varying vaccine duration
#' @param model_out Model object
#' @export
get_model_likelihood.vacc_durR_nimue_simulation <- function(model_out){
  excess_log_likelihood_vaccine
}
#' Method to calculate log-likelihood for each model type for a nimue model for LMIC work
#' @param model_out Model object
#' @export
get_model_likelihood.lmic_nimue_simulation <- function(model_out){
  calc_loglikelihood_variant
}
#' Method to calculate log-likelihood for each model type for a nimue model for LMIC work with boosters
#' @param model_out Model object
#' @export
get_model_likelihood.lmic_booster_nimue_simulation <- function(model_out){
  calc_loglikelihood_booster
}
