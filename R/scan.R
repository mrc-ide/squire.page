#' Fit a nimue model to the given data using a scan fitting method
#'
#' Taken from OJs code. Iteratively adjusts the Rt_rw parameters based on
#' maximising likelihood of a short period of time in the data.
#'
#' Only works with excess fits for now.
#'
#' Uses furrr for simulations later so can be run in parallel with
#' future::plan(future::multisession(), .cleanup = TRUE)
#'
#' @param res A nimue simulation object
#' @param data The data to fit to, for standard fits must have date and death
#' columns; for lmic fits must have date, deaths and cases; and for excess
#' mortality it must have deaths, week_start and week_end. Default = NULL, which uses the data in the fit
#' @param width The algorithm will check values up to +- width of the current
#' value
#' @param n_span How many values between the width to check
#' @param width_end The final width to check, so that the algorithm to move to
#' a finer resolution as the proceeds, default NULL means the same as width
#' @param repeats How many times to iterate through the process
#' @return Model object tuned to the data
#' @export
# scan_fit <- function(res, data = NULL, width = 2, n_span = 8, width_end = NULL, repeats = 5) {
#   if(is.null(width_end)){
#     widths = rep(width, repeats)
#   } else {
#     widths = round(seq(width, width_end, length.out = repeats))
#   }
#   if(is.null(data)){
#     data <- res$pmcmc_results$inputs$data
#   }
#   #calculate the data rows, from time kicks in to plus 20 days
#   if("excess_nimue_simulation" %in% class(res)){
#     dates_Rt <- res$interventions$date_Rt_change
#     #now get the relevant dates
#     data_indexes <- suppressWarnings(
#       lapply(dates_Rt, function(date){
#         start_index <- min(which(data$week_start > date))
#         end_index <- max(
#           max(which(data$week_start < date + 21)),
#           1)
#         seq(start_index, end_index, 1)
#       })
#     )
#   } else {
#     stop("This function is currently step up to work excess fits only")
#   }
#   #get the likelihood function
#   likelihood_func <- get_model_likelihood(res)
#   #get pars from the replicates for now
#   pars <- res$replicate_parameters
#   #run through repeats and parameters
#   for(rep in seq.int(repeats)){
#     message(paste0("Repeat: ", rep))
#     # loop over each parameter and scan values
#     for(i in seq_along(data_indexes)) {
#       message(paste0("Parameter:", i))
#       data_i <- data[data_indexes[[i]],]
#       name <- paste0("Rt_rw_", i)
#       change <- c(
#         setdiff(seq(-widths[rep], widths[rep], length.out = n_span), 0),
#         0
#         )
#       #calculate likelihood
#       ll_i <- furrr::future_map_dbl(
#         .x = change,
#         .f = function(change){
#           pars_temp <- pars
#           pars_temp[, name] <- pars[, name] + change
#           #pars[, name]  <- pars[, name] + change
#           purrr::map_dbl(
#             .x = seq.int(nrow(pars)),
#             .f = function(x){
#               likelihood_func(
#                 pars = pars_temp[x,],
#                 data = data_i,
#                 squire_model = res$pmcmc_results$inputs$squire_model,
#                 model_params = res$pmcmc_results$inputs$model_params,
#                 pars_obs = res$pmcmc_results$inputs$pars_obs,
#                 n_particles = 2, forecast_days = 0, return = "ll",
#                 Rt_args = res$pmcmc_results$inputs$Rt_args,
#                 interventions = res$pmcmc_results$inputs$interventions
#               )$log_likelihood
#             }
#           ) %>%
#             mean()
#         },
#         .options = furrr::furrr_options(seed = TRUE)
#         #this requires more thought, since now might not be reproducible,
#         #however the effect of the noise is so small that it has almost
#         #no impact and so should be fine for now
#       )
#       #update parameter
#       pars[, name] <- pars[, name] + change[which.max(ll_i)]
#     }
#   }
#   #now we need to update our replicates with the new pars and an output
#   res$replicate_parameters <- pars
#   #make the pars into a pars_list and then simulate and output
#   #NOTE: as of now this leaves the chain itself as it is
#   res <- generate_draws(res, pars.list = NULL, draws = NULL)
#   return(res)
# }
scan_fit <- function(res, data = NULL, width = 2, n_span = 8, width_end = NULL, repeats = 5) {
  if(is.null(width_end)){
    widths = rep(width, repeats)
  } else {
    widths = seq(width, width_end, length.out = repeats)
  }
  if(is.null(data)){
    data <- res$pmcmc_results$inputs$data
  }
  #calculate the data rows, from time kicks in to plus 20 days
  if("excess_nimue_simulation" %in% class(res)){
    dates_Rt <- res$interventions$date_Rt_change
    #now get the relevant dates
    data_indexes <- suppressWarnings(
      lapply(dates_Rt, function(date){
        start_index <- max(max(which(data$week_start < date)), 1)
        end_index <- max(
          max(which(data$week_start <= date + 21)),
          1)
        seq(start_index, end_index, 1)
      })
    )
  } else {
    stop("This function is currently step up to work excess fits only")
  }
  #get the likelihood function
  likelihood_func <- get_model_likelihood(res)
  #get the average pars from the replicates
  pars_names <- setdiff(names(res$replicate_parameters), c("start_date", "R0"))
  avg_pars <- colMeans(res$replicate_parameters[, pars_names])
  #convert to data frame
  pars <- as.data.frame(t(avg_pars))
  #add start date and R0
  pars$start_date <- stats::median(res$replicate_parameters[,"start_date"])
  pars$R0 <- mean(res$replicate_parameters[,"R0"])
  #run through repeats and parameters
  for(rep in seq.int(repeats)){
    message(paste0("Repeat: ", rep))
    # loop over each parameter and scan values
    for(i in seq_along(data_indexes)) {
      if(length(data_indexes[[i]]) > 2){

        #message(paste0("Parameter:", i))
        data_i <- data[data_indexes[[i]],]
        name <- paste0("Rt_rw_", i)
        change <- c(
          setdiff(seq(-widths[rep], widths[rep], length.out = n_span), 0),
          0
        )
        #calculate likelihood
        ll_i <- purrr::map_dbl(
          .x = change,
          .f = function(change){
            pars_temp <- pars
            pars_temp[1, name] <- pars[1, name] + change
            likelihood_func(
              pars = pars_temp[1, ],
              data = data_i,
              squire_model = res$pmcmc_results$inputs$squire_model,
              model_params = res$pmcmc_results$inputs$model_params,
              pars_obs = res$pmcmc_results$inputs$pars_obs,
              n_particles = 2, forecast_days = 0, return = "ll",
              Rt_args = res$pmcmc_results$inputs$Rt_args,
              interventions = res$pmcmc_results$inputs$interventions
            )$log_likelihood +
              res$pmcmc_results$inputs$prior(pars_temp)
          }
        )
        #update parameter
        pars[1, name] <- pars[1, name] + change[which.max(ll_i)]
      }
    }
  }
  #now we need to update our replicates with the new pars and an output
  #scale replicates to get their mean to our current replicate
  #calculate difference from average
  avg_diff <- avg_pars[pars_names] - unlist(pars[1, pars_names])[pars_names]
  #scale just the replicate for now
  res$replicate_parameters[, pars_names] <- sweep(
    res$replicate_parameters[, pars_names],
    2,
    avg_diff[pars_names],
    FUN = "-"
  )
  #make the pars into a pars_list and then simulate and output
  res <- generate_draws(res, pars.list = NULL, draws = NULL)
  #this does not effect the chain itself so we scale that now and recalculate
  #posterior
  res$pmcmc_results$results[, pars_names] <-
    sweep(
      res$pmcmc_results$results[, pars_names],
      2,
      avg_diff[pars_names],
      FUN = "-"
    )
  #we do not update the posterior/priors/likelihood for now as this would take too long
  # #update prior, likelihood and posterior
  # for(row in 1:nrow(res$pmcmc_results$results)) {
  #   row_pars <- res$pmcmc_results$results[row, ]
  #   prior <- res$pmcmc_results$inputs$prior(
  #     row_pars
  #   )
  #   #make start date a date
  #   row_pars$start_date <- get_data_start_date(res) + row_pars$start_date
  #   likelihood <- likelihood_func(
  #     pars = row_pars,
  #     data = data,
  #     squire_model = res$pmcmc_results$inputs$squire_model,
  #     model_params = res$pmcmc_results$inputs$model_params,
  #     pars_obs = res$pmcmc_results$inputs$pars_obs,
  #     n_particles = 2, forecast_days = 0, return = "ll",
  #     Rt_args = res$pmcmc_results$inputs$Rt_args,
  #     interventions = res$pmcmc_results$inputs$interventions
  #   )$log_likelihood
  #   res$pmcmc_results$results[row,
  #                             c("log_prior",
  #                               "log_likelihood",
  #                               "log_posterior")] <- c(
  #                                 prior, likelihood, prior + likelihood
  #                               )
  # }
  return(res)
}
