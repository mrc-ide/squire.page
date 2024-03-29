get_model_likelihood <- function(model_out){
  UseMethod("get_model_likelihood")
}
get_model_likelihood.nimue_simulation <- function(model_out){
  squire:::calc_loglikelihood
}
get_model_likelihood.excess_nimue_simulation <- function(model_out){
  excess_log_likelihood
}
get_model_likelihood.vacc_durR_nimue_simulation <- function(model_out){
  excess_log_likelihood_vaccine
}
get_model_likelihood.lmic_nimue_simulation <- function(model_out){
  calc_loglikelihood_variant
}
get_model_likelihood.lmic_booster_nimue_simulation <- function(model_out){
  calc_loglikelihood_booster
}
