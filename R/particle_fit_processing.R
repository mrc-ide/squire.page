#' Trim an Rt Optimised output
#'
#' Removes poorly fitted trajectories based upon a predefined cutoff.
#' For each optimised trajectory this calculates the least squares error from
#' the observed data, then scales those errors so that the error for estimating
#' 0's is 1 then rejects any trajectories with an error higher than the p_error
#' value.
#'
#'
#' @param out An object of type rt_optimised.
#' @param p_error Cutoff value, see description.
#'
#' @return An object of type rt_optimised_trimmed, (model type).
#'
#' @export
trim_rt_optimise <- function(out, p_error){
  #checks
  if(is.null(out$output)){
    stop("out must have generated model outputs, use generate_draws to get these values.")
  }
  #get data deaths
  data_deaths <- out$inputs$data$deaths
  t_start <- out$inputs$data$t_start
  t_end <- out$inputs$data$t_end
  err_func <- function(model_deaths, data_deaths){
    sqrt(sum((model_deaths - data_deaths)^2))
  }
  #for each sample
  errors <- get_cumulative_deaths_rt_optimise_output(out) %>%
    #is it within the error
    purrr::map_dbl(function(cumulative_deaths){
      deaths <- cumulative_deaths[t_end] - cumulative_deaths[t_start]
      err_func(deaths, data_deaths)
    })
  #rescale so that the error of all 0's is 1
  error_0 <- err_func(rep(0, length(data_deaths)), data_deaths)
  errors <- errors/error_0
  #if errors are below the given p_error value we keep them
  good_fits_index <- which(errors < p_error)
  bad_fits_index <- setdiff(seq_along(out$samples), good_fits_index)
  #split samples up, if any don't pass
  if(length(bad_fits_index) > 0){
    out$excluded$samples <- out$samples[bad_fits_index]
    out$excluded$output <- out$output[, , bad_fits_index]
    #ensure its a 3d array
    dim(out$excluded$output) <- c(dim(out$excluded$output)[1:2], length(bad_fits_index))
    out$samples <- out$samples[good_fits_index]
    out$output <- out$output[, , good_fits_index]
    class(out) <- c("rt_optimised_trimmed", class(out))
  }
  out
}

#' Function to get the likelihoods of attached odin model outputs
#' @noRd
calc_likelihoods_rt_optimise_output <- function(out){
  #get data deaths
  data_deaths <- out$inputs$data$deaths
  t_start <- out$inputs$data$t_start
  t_end <- out$inputs$data$t_end
  #for each sample
  get_cumulative_deaths_rt_optimise_output(out) %>%
    purrr::map_dbl(function(cumulative_deaths){
      deaths <- cumulative_deaths[t_end] - cumulative_deaths[t_start]
      #get the log likelihood
      ll_negative_binomial(deaths, data_deaths, out$inputs$k)
    })
}

#' Function to get the deaths from attached odin model outputs
#' @noRd
get_cumulative_deaths_rt_optimise_output <- function(out){
  #for each sample
  purrr::map(seq_len(dim(out$output)[3]), function(slice){
    #extract deaths
    cumulative_deaths <- out$output[, stringr::str_detect(colnames(out$output[, , slice]), "D\\["), slice] %>%
      rowSums()
  })
}
