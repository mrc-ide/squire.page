get_data_end_date <- function(model_out){
  UseMethod("get_data_end_date")
}
get_data_end_date.nimue_simulation <- function(model_out){
  max(model_out$pmcmc_results$inputs$data$date)
}
get_data_end_date.excess_nimue_simulation <- function(model_out){
  max(model_out$pmcmc_results$inputs$data$week_end)
}

get_data_end_date_inner <- function(model_out){
  UseMethod("get_data_end_date_inner")
}
get_data_end_date_inner.nimue_simulation <- function(model_out){
  max(model_out$pmcmc_results$inputs$data$date)
}
get_data_end_date_inner.excess_nimue_simulation <- function(model_out){
  max(model_out$pmcmc_results$inputs$data$week_start)
}

get_data_start_date <- function(model_out){
  UseMethod("get_data_start_date")
}
get_data_start_date.nimue_simulation <- function(model_out){
  min(model_out$pmcmc_results$inputs$data$date)
}
get_data_start_date.excess_nimue_simulation <- function(model_out){
  min(model_out$pmcmc_results$inputs$data$week_start)
}

get_dates <- function(model_out){
  UseMethod("get_dates")
}
get_dates.nimue_simulation <- function(model_out){
  model_out$pmcmc_results$inputs$data$date
}
get_dates.excess_nimue_simulation <- function(model_out){
  model_out$pmcmc_results$inputs$data$week_start
}
get_dates_greater <- function(model_out){
  UseMethod("get_dates_greater")
}
get_dates_greater.nimue_simulation <- function(model_out){
  model_out$pmcmc_results$inputs$data$date
}
get_dates_greater.excess_nimue_simulation <- function(model_out){
  model_out$pmcmc_results$inputs$data$week_end
}
