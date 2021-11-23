#' Scale a Nimue MCMC output by the number of reported cases
#'
#' Performs a simple grid search, adding or subtracting from the final fitted
#' effect (in all replicates of the model object). Chooses the value that
#' minimises the error from the total number of infections in the final
#' \code{grad_dur} days of the epidemic, adjusted by a reporting fraction
#' calculated over the previous few days (14 by default, determined by
#' \code{Rt_rw_duration}). A legacy function and is no longer used.
#'
#' @param out A nimue model object.
#' @param n_particles How many replicates to draw from the adjust object.
#' @param grad_dur How many days into the past to adjust.
#' @return A nimue model object, with the adjusted replicate parameters.
#' @export
generate_draws_pmcmc_nimue_case_fitted <- function(out, n_particles = 10, grad_dur = 21) {

  pmcmc <- out$pmcmc_results
  n_chains <- max(length(out$pmcmc_results$chains), 1)
  burnin <- round(out$pmcmc_results$inputs$n_mcmc/10)
  squire_model <- out$pmcmc_results$inputs$squire_model
  replicates <- dim(out$output)[3]
  forecast <- 0
  country <- out$parameters$country
  population <- out$parameters$population
  interventions <- out$interventions
  data <- out$pmcmc_results$inputs$data
  rw_dur <- out$pmcmc_results$inputs$Rt_args$Rt_rw_duration

  #--------------------------------------------------------
  # Section 1 # what is our predicted gradient
  #--------------------------------------------------------

  # first what is the model predicted infections
  infections <- nimue_format(out, "infections", date_0 = max(data$date))
  infections_end <- infections %>% dplyr::filter(date > (max(data$date) - grad_dur) & date <= (max(data$date))) %>%
    dplyr::group_by(date) %>% dplyr::summarise(y = stats::median(.data$y))

  infections_pre_end <- infections %>%
    dplyr::filter(.data$date > (max(data$date) - grad_dur - rw_dur) & date <= (max(data$date) - grad_dur) ) %>%
    dplyr::group_by(.data$date) %>% dplyr::summarise(y = stats::median(.data$y))

  # and the observed cases
  cases_end <- utils::tail(data$cases, grad_dur)
  cases_pre_end <- utils::head(utils::tail(data$cases, grad_dur+rw_dur), rw_dur)

  # get these gradients
  get_infs <- function(x) {
    sum(x, na.rm = TRUE)
  }

  pred_infs_end <- get_infs(infections_end$y)
  pred_infs_pre_end <- get_infs(infections_pre_end$y)

  des_infs_end <- get_infs(cases_end)
  des_infs_pre_end <- get_infs(cases_pre_end)

  # if there are less than 100 cases in both windowns then don't bother
  if(des_infs_end > 100 && des_infs_pre_end > 100) {

    ca_infs_frac <-  pred_infs_pre_end / des_infs_pre_end

    # desired model predictd final infs
    wanted_infs <- des_infs_end * ca_infs_frac

    # if actual infs available
    if(!is.nan(wanted_infs) || !is.na(wanted_infs) || !is.infinite(wanted_infs)) {

      # do we need to go up or down
      if(wanted_infs < pred_infs_end) {
        alters <- seq(0.025, 0.2, 0.025)
      } else {
        alters <- seq(-0.025, -0.15, -0.025) # more conservative on the way up
      }

      # store our grads
      ans <- alters
      if ("chains" %in% names(out$pmcmc_results)) {
        last_rw <- ncol(out$pmcmc_results$chains$chain1$results) - 3
      } else {
        last_rw <- ncol(out$pmcmc_results$results) - 3
      }

      #--------------------------------------------------------
      # Section 2 # # find best grad correction
      #--------------------------------------------------------

      for(alt in seq_along(alters)) {

        message(alt)

        if ("chains" %in% names(out$pmcmc_results)) {
          for(ch in seq_along(out$pmcmc_results$chains)) {
            out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] + alters[alt]
          }
        } else {
          out$pmcmc_results$results[,last_rw] <- out$pmcmc_results$results[,last_rw] + alters[alt]
        }

        pmcmc_samples <- squire:::sample_pmcmc(pmcmc_results = out$pmcmc_results,
                                               burnin = burnin,
                                               n_chains = n_chains,
                                               n_trajectories = replicates,
                                               n_particles = n_particles,
                                               forecast_days = forecast)

        dimnms <- dimnames(pmcmc_samples$trajectories)

        # then let's create the output that we are going to use
        names(pmcmc_samples)[names(pmcmc_samples) == "trajectories"] <- "output"
        dimnames(pmcmc_samples$output) <- list(dimnames(pmcmc_samples$output)[[1]], dimnames(out$output)[[2]], NULL)
        out$output <- pmcmc_samples$output

        # and adjust the time as before
        full_row <- match(0, apply(out$output[,"time",],2,function(x) { sum(is.na(x)) }))
        saved_full <- out$output[,"time",full_row]
        for(i in seq_len(replicates)) {
          na_pos <- which(is.na(out$output[,"time",i]))
          full_to_place <- saved_full - which(rownames(out$output) == as.Date(max(data$date))) + 1L
          if(length(na_pos) > 0) {
            full_to_place[na_pos] <- NA
          }
          out$output[,"time",i] <- full_to_place
        }

        infections <- nimue_format(out, "infections", date_0 = max(data$date))
        this_infs <- infections %>% dplyr::filter(date > (max(data$date) - grad_dur) & date <= (max(data$date))) %>%
          dplyr::group_by(date) %>% dplyr::summarise(y = stats::median(.data$y))

        ans[alt] <- get_infs(this_infs$y)

        # put our chains back to normal
        if ("chains" %in% names(out$pmcmc_results)) {
          for(ch in seq_along(out$pmcmc_results$chains)) {
            out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] - alters[alt]
          }
        } else {
          out$pmcmc_results$results[,last_rw] <- out$pmcmc_results$results[,last_rw] - alters[alt]
        }

      }


      # adapt our whole last chain accordingly
      alts <- which.min(abs(ans-wanted_infs))
      if ("chains" %in% names(out$pmcmc_results)) {
        for(ch in seq_along(out$pmcmc_results$chains)) {
          out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] + alters[alts]
        }
      } else {
        out$pmcmc_results$results[,last_rw] <- out$pmcmc_results$results[,last_rw] + alters[alt]
      }

    }

  }

  #--------------------------------------------------------
  # Section 3 of pMCMC Wrapper: Sample PMCMC Results
  #--------------------------------------------------------
  pmcmc_samples <- squire:::sample_pmcmc(pmcmc_results = out$pmcmc_results,
                                         burnin = burnin,
                                         n_chains = n_chains,
                                         n_trajectories = replicates,
                                         n_particles = n_particles,
                                         forecast_days = forecast)

  #--------------------------------------------------------
  # Section 4 of pMCMC Wrapper: Tidy Output
  #--------------------------------------------------------
  dimnms <- dimnames(pmcmc_samples$trajectories)

  # then let's create the output that we are going to use
  names(pmcmc_samples)[names(pmcmc_samples) == "trajectories"] <- "output"
  dimnames(pmcmc_samples$output) <- list(dimnames(pmcmc_samples$output)[[1]], dimnames(out$output)[[2]], NULL)
  out$output <- pmcmc_samples$output
  out$replicate_parameters <- pmcmc_samples$sampled_PMCMC_Results

  # and adjust the time as before
  full_row <- match(0, apply(out$output[,"time",],2,function(x) { sum(is.na(x)) }))
  saved_full <- out$output[,"time",full_row]
  for(i in seq_len(replicates)) {
    na_pos <- which(is.na(out$output[,"time",i]))
    full_to_place <- saved_full - which(rownames(out$output) == as.Date(max(data$date))) + 1L
    if(length(na_pos) > 0) {
      full_to_place[na_pos] <- NA
    }
    out$output[,"time",i] <- full_to_place
  }

  return(out)

}

#' Scale a Squire MCMC output by the number of reported cases
#'
#' Performs a simple grid search, adding or subtracting from the final fitted
#' effect (in all replicates of the model object). Chooses the value that
#' minimises the error from the total number of infections in the final
#' \code{grad_dur} days of the epidemic, adjusted by a reporting fraction
#' calculated over the previous few days (14 by default, determined by
#' \code{Rt_rw_duration}). A legacy function and is no longer used.
#'
#' @param out A squire model object.
#' @param n_particles How many replicates to draw from the adjust object.
#' @param grad_dur How many days into the past to adjust.
#' @return A squire model object, with the adjusted replicate parameters.
#' @export
#'@export
generate_draws_pmcmc_case_fitted <- function(out, n_particles = 10, grad_dur = 21) {

  pmcmc <- out$pmcmc_results
  n_chains <- length(out$pmcmc_results$chains)
  burnin <- round(out$pmcmc_results$inputs$n_mcmc/10)
  squire_model <- out$pmcmc_results$inputs$squire_model
  replicates <- dim(out$output)[3]
  forecast <- 0
  country <- out$parameters$country
  population <- out$parameters$population
  interventions <- out$interventions
  data <- out$pmcmc_results$inputs$data
  rw_dur <- out$pmcmc_results$inputs$Rt_args$Rt_rw_duration

  #--------------------------------------------------------
  # Section 1 # what is our predicted gradient
  #--------------------------------------------------------

  # first what is the model predicted infections
  infections <- squire::format_output(out, "infections", date_0 = max(data$date))
  infections_end <- infections %>% dplyr::filter(date > (max(data$date) - grad_dur) & date <= (max(data$date))) %>%
    dplyr::group_by(date) %>% dplyr::summarise(y = stats::median(.data$y))

  infections_pre_end <- infections %>%
    dplyr::filter(date > (max(data$date) - grad_dur - rw_dur) & date <= (max(data$date) - grad_dur) ) %>%
    dplyr::group_by(date) %>% dplyr::summarise(y = stats::median(.data$y))

  # and the observed cases
  cases_end <- utils::tail(data$cases, grad_dur)
  cases_pre_end <- utils::head(utils::tail(data$cases, grad_dur+rw_dur), rw_dur)

  # get these gradients
  get_infs <- function(x) {
    sum(x, na.rm = TRUE)
  }

  pred_infs_end <- get_infs(infections_end$y)
  pred_infs_pre_end <- get_infs(infections_pre_end$y)

  des_infs_end <- get_infs(cases_end)
  des_infs_pre_end <- get_infs(cases_pre_end)

  # if there are less than 100 cases in both windowns then don't bother
  if(des_infs_end > 100 && des_infs_pre_end > 100) {

    ca_infs_frac <-  pred_infs_pre_end / des_infs_pre_end

    # desired model predictd final infs
    wanted_infs <- des_infs_end * ca_infs_frac

    # if actual infs available
    if(!is.nan(wanted_infs) || !is.na(wanted_infs) || !is.infinite(wanted_infs)) {

      index <- squire:::odin_index(out$model)
      index$n_E2_I <- seq(utils::tail(unlist(index),1)+1, utils::tail(unlist(index),1)+length(index$S),1)
      index$delta_D <- seq(utils::tail(unlist(index),1)+1, utils::tail(unlist(index),1)+length(index$S),1)

      # do we need to go up or down
      if(wanted_infs < pred_infs_end) {
        alters <- seq(0.025, 0.175, 0.025)
      } else {
        alters <- seq(-0.025, -0.125, -0.025) # more conservative on the way up
      }

      # store our grads
      ans <- alters
      last_rw <- ncol(out$pmcmc_results$chains$chain1$results) - 3

      # for later
      all_case_compartments <- unlist(
        index[c("IMild", "ICase1", "ICase2", "IOxGetLive1", "IOxGetLive2",
                "IOxGetDie1", "IOxGetDie2", "IOxNotGetLive1", "IOxNotGetLive2",
                "IOxNotGetDie1", "IOxNotGetDie2", "IMVGetLive1", "IMVGetLive2",
                "IMVGetDie1", "IMVGetDie2", "IMVNotGetLive1", "IMVNotGetLive2",
                "IMVNotGetDie1", "IMVNotGetDie2", "IRec1", "IRec2", "R", "D")])

      #--------------------------------------------------------
      # Section 2 # # find best grad correction
      #--------------------------------------------------------

      for(alt in seq_along(alters)) {

        message(alt)

        for(ch in seq_along(out$pmcmc_results$chains)) {
          out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] + alters[alt]
        }

        pmcmc_samples <- squire:::sample_pmcmc(pmcmc_results = out$pmcmc_results,
                                               burnin = burnin,
                                               n_chains = n_chains,
                                               n_trajectories = replicates,
                                               n_particles = n_particles,
                                               forecast_days = forecast)

        dimnms <- dimnames(pmcmc_samples$trajectories)

        # make e2_i space
        dimnms[[2]] <- c(dimnms[[2]], paste0("n_E2_I[", seq_len(length(index$S)),"]"))
        new_data <- array(data = 0,
                          dim = c(dim(pmcmc_samples$trajectories) + c(0, length(index$n_E2_I), 0)),
                          dimnames = dimnms)

        new_data[, seq_len(dim(pmcmc_samples$trajectories)[2]), ] <- pmcmc_samples$trajectories
        pmcmc_samples$trajectories <- new_data

        # make D space
        dimnms <- dimnames(pmcmc_samples$trajectories)
        dimnms[[2]] <- c(dimnms[[2]], paste0("delta_D[", seq_len(length(index$S)),"]"))
        new_data <- array(data = 0,
                          dim = c(dim(pmcmc_samples$trajectories) + c(0, length(index$delta_D), 0)),
                          dimnames = dimnms)

        new_data[, seq_len(dim(pmcmc_samples$trajectories)[2]), ] <- pmcmc_samples$trajectories
        pmcmc_samples$trajectories <- new_data
        nt <- nrow(pmcmc_samples$trajectories)

        # are the steps not 1 apart? if so we need to sum the incident variables (infecions/deaths)
        if (out$parameters$day_return || !squire:::odin_is_discrete(out$model)) {

          # assign the infections
          for(i in seq_along(out$parameters$population)) {
            collect <- vapply(1:out$parameters$replicates, function(j) {
              pos <- seq(i, length(index$cum_infs), by = length(out$parameters$population))
              pos <- index$cum_infs[pos]
              diff(pmcmc_samples$trajectories[,pos,j])
            }, FUN.VALUE = numeric(nt-1))
            pmcmc_samples$trajectories[1+seq_len(nt-1),index$n_E2_I[i],] <- collect
          }

          # assign the deaths
          for(i in seq_along(out$parameters$population)) {
            collect <- vapply(1:out$parameters$replicates, function(j) {
              pos <- seq(i, length(index$D), by = length(out$parameters$population))
              pos <- index$D[pos]
              diff(pmcmc_samples$trajectories[,pos,j])
            }, FUN.VALUE = numeric(nt-1))
            pmcmc_samples$trajectories[1+seq_len(nt-1),index$delta_D[i],] <- collect
          }

        }

        this_infs <- as.numeric(rowMeans(utils::tail(matrix(unlist(lapply(seq_len(replicates), function(i) {
          rowSums(pmcmc_samples$trajectories[,index$n_E2_I,i])
        })), ncol = replicates),grad_dur)))

        ans[alt] <- get_infs(this_infs)


        for(ch in seq_along(out$pmcmc_results$chains)) {
          out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] - alters[alt]
        }

      }


      # adapt our whole last chain accordingly
      alts <- which.min(abs(ans-wanted_infs))
      for(ch in seq_along(out$pmcmc_results$chains)) {
        out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] + alters[alts]
      }

    }

  }

  # set up now to do the stochastic draws
  out$pmcmc_results$inputs$squire_model <- squire::explicit_model()
  out$pmcmc_results$inputs$model_params$dt <- 0.02

  #--------------------------------------------------------
  # Section 3 of pMCMC Wrapper: Sample PMCMC Results
  #--------------------------------------------------------
  pmcmc_samples <- squire:::sample_pmcmc(pmcmc_results = out$pmcmc_results,
                                         burnin = burnin,
                                         n_chains = n_chains,
                                         n_trajectories = replicates,
                                         n_particles = n_particles,
                                         forecast_days = forecast)

  #--------------------------------------------------------
  # Section 4 of pMCMC Wrapper: Tidy Output
  #--------------------------------------------------------

  # create a fake run object and fill in the required elements
  r <- out$pmcmc_results$inputs$squire_model$run_func(country = out$parameters$country,
                                                      contact_matrix_set = out$pmcmc_results$inputs$model_params$contact_matrix_set,
                                                      tt_contact_matrix = out$pmcmc_results$inputs$model_params$tt_matrix,
                                                      hosp_bed_capacity = out$pmcmc_results$inputs$model_params$hosp_bed_capacity,
                                                      tt_hosp_beds = out$pmcmc_results$inputs$model_params$tt_hosp_beds,
                                                      ICU_bed_capacity = out$pmcmc_results$inputs$model_params$ICU_bed_capacity,
                                                      tt_ICU_beds = out$pmcmc_results$inputs$model_params$tt_ICU_beds,
                                                      population = out$pmcmc_results$inputs$population,
                                                      replicates = 1,
                                                      day_return = TRUE,
                                                      time_period = nrow(pmcmc_samples$trajectories),
                                                      dur_R = out$pmcmc_results$inputs$model_params$dur_R)

  # and add the parameters that changed between each simulation, i.e. posterior draws
  r$replicate_parameters <- pmcmc_samples$sampled_PMCMC_Results

  # as well as adding the pmcmc chains so it's easy to draw from the chains again in the future
  r$pmcmc_results <- out$pmcmc_results

  # then let's create the output that we are going to use
  names(pmcmc_samples)[names(pmcmc_samples) == "trajectories"] <- "output"
  dimnames(pmcmc_samples$output) <- list(dimnames(pmcmc_samples$output)[[1]], dimnames(r$output)[[2]], NULL)
  r$output <- pmcmc_samples$output

  # and adjust the time as before
  full_row <- match(0, apply(r$output[,"time",],2,function(x) { sum(is.na(x)) }))
  saved_full <- r$output[,"time",full_row]
  for(i in seq_len(replicates)) {
    na_pos <- which(is.na(r$output[,"time",i]))
    full_to_place <- saved_full - which(rownames(r$output) == as.Date(max(data$date))) + 1L
    if(length(na_pos) > 0) {
      full_to_place[na_pos] <- NA
    }
    r$output[,"time",i] <- full_to_place
  }

  # second let's recreate the output
  r$model <- pmcmc_samples$inputs$squire_model$odin_model(
    user = pmcmc_samples$inputs$model_params, unused_user_action = "ignore"
  )

  # we will add the interventions here so that we know what times are needed for projection
  r$interventions <- interventions

  # and fix the replicates
  r$parameters$replicates <- replicates
  r$parameters$time_period <- as.numeric(diff(as.Date(range(rownames(r$output)))))
  r$parameters$dt <- out$pmcmc_results$inputs$model_params$dt

  return(r)

}
