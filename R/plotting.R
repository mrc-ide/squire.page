#' A diagnostic plot for looking at the infections-based fitting
#' @param out model object
#'@export
compare_adjustment_plot <- function(out){
  #get the model infections
  plotting_data <- nimue_format(out, var_select = "infections") %>%
    dplyr::filter(t >= (max(.data$t, na.rm = TRUE) - out$pmcmc_results$inputs$pars_obs$cases_days - out$pmcmc_results$inputs$pars_obs$cases_reporting)) %>%
    dplyr::rename(model_infections = "y") %>%
    dplyr::select(replicate, .data$t, .data$model_infections) %>%
    dplyr::mutate(
      reporting = .data$t < max(.data$t, na.rm = TRUE) - out$pmcmc_results$inputs$pars_obs$cases_days
    ) %>%
    #add the real data
    dplyr::left_join(
      get_data(out) %>%
        dplyr::mutate(t = as.numeric(.data$date - max(.data$date))) %>%
        dplyr::select(!.data$deaths),
      by = "t"
    ) %>%
    #calculate reporting fraction
    dplyr::group_by(replicate) %>%
    dplyr::mutate(
      reporting_fraction = mean(
        dplyr::if_else(.data$reporting,
                       .data$cases/.data$model_infections,
                       as.numeric(NA)),
        na.rm = TRUE
      )
    ) %>%
    #convert cases to the ones we fit too
    dplyr::mutate(
      adjusted_cases = .data$cases/.data$reporting_fraction
    )
  #plot
  ggplot2::ggplot(data = plotting_data, ggplot2::aes(x = .data$date, y = .data$adjusted_cases)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = max(.data$date)-out$pmcmc_results$inputs$pars_obs$cases_days),
               linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.25) +
    ggplot2::geom_line(
      inherit.aes = FALSE,
      data = plotting_data %>%
        dplyr::group_by(.data$date) %>%
        dplyr::summarise(
          infections_med = stats::median(.data$model_infections , na.rm = TRUE)
        ),
      ggplot2::aes(x = .data$date, y = .data$infections_med),
      colour = "red"
    ) +
    ggplot2::geom_ribbon(
      inherit.aes = FALSE,
      data = plotting_data %>%
        dplyr::group_by(.data$date) %>%
        dplyr::summarise(
          infections_025 = stats::quantile(.data$model_infections , 0.025, na.rm = TRUE),
          infections_975 = stats::quantile(.data$model_infections , 0.975, na.rm = TRUE)
        ),
      ggplot2::aes(x = .data$date, ymin = .data$infections_025,
          ymax = .data$infections_975),
      alpha =0.25,
      fill = "red"
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Date", y = "Infections, cases adjusted\nby reporting fraction") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
                   axis.title.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black")
    ) + ggplot2::theme(legend.position = "top",
                       legend.justification = c(0,1),
                       legend.direction = "horizontal") +
    ggplot2::lims(y = c(0, max(c(plotting_data$adjusted_cases,
                        plotting_data$model_infections))))
}
#' An S3 generic to calculate the proportion of the population immune
#' @param out squire/nimue model object
#' @param max_date date to generate immunity ratio upto, default = NULL gets the max_date from out's data
#' @param vaccine should the function account for vaccinations, default = FALSE
#' @export
get_immunity_ratios <- function(out, max_date = NULL, vaccine = FALSE){
  UseMethod("get_immunity_ratios")
}
#' An S3 method to calculate the proportion of the population immune
#' @inheritParams get_immunity_ratios
#' @export
get_immunity_ratios.default <- function(out, max_date = NULL, vaccine = FALSE) {
  if(is.null(out$pmcmc_results)){
    mixing_matrix <- squire:::process_contact_matrix_scaled_age(
      out$odin_parameters$contact_matrix_set[[1]],
      out$odin_parameters$population
    )

    dur_ICase <- out$parameters$dur_ICase
    dur_IMild <- out$parameters$dur_IMild
    prob_hosp <- out$parameters$prob_hosp

    # assertions
    squire:::assert_single_pos(dur_ICase, zero_allowed = FALSE)
    squire:::assert_single_pos(dur_IMild, zero_allowed = FALSE)
    squire:::assert_numeric(prob_hosp)
    squire:::assert_numeric(mixing_matrix)
    squire:::assert_square_matrix(mixing_matrix)
    squire:::assert_same_length(mixing_matrix[,1], prob_hosp)

    if(sum(is.na(prob_hosp)) > 0) {
      stop("prob_hosp must not contain NAs")
    }

    if(sum(is.na(mixing_matrix)) > 0) {
      stop("mixing_matrix must not contain NAs")
    }

    index <- squire:::odin_index(out$model)
    pop <- out$parameters$population

    if(is.null(max_date)) {
      max_date <- max(get_dates(out))
    }
    t_now <- which(as.Date(rownames(out$output)) == max_date)

    if(vaccine){
      # make the vaccine time changing args
      nrs <- t_now
      vei <- lapply(seq_len(nrow(out$odin_parameters$vaccine_efficacy_infection)),
                    function(x) {out$odin_parameters$vaccine_efficacy_infection[x,,]}
      )
      phl <- lapply(seq_len(nrow(out$odin_parameters$prob_hosp)),
                    function(x) {out$odin_parameters$prob_hosp[x,,]}
      )

      if(nrs > length(vei)) {
        vei_full <- c(rep(list(vei[[1]]),nrs - length(vei)), vei)
        phl_full <- c(rep(list(phl[[1]]),nrs - length(phl)), phl)
      } else {
        vei_full <- utils::tail(vei, nrs)
        phl_full <- utils::tail(phl, nrs)
      }

      # prop susceptible
      prop_susc <- lapply(seq_len(dim(out$output)[3]), function(x) {

        susceptible <- array(
          out$output[seq_len(t_now),index$S,x],
          dim=c(t_now, dim(index$S))
        )

        # We divide by the total population
        prop_susc <- sweep(susceptible, 2, pop, FUN='/')

        # We multiply by the effect of vaccines on onward infectiousness
        prop_susc <- vapply(
          seq_len(nrow(prop_susc)),
          FUN = function(i){prop_susc[i,,]*vei_full[[i]]},
          FUN.VALUE = prop_susc[1,,]
        )

        prop_susc <- aperm(prop_susc, c(3,1,2))

        return(prop_susc)
      } )

      relative_R0_by_age <- lapply(phl_full, function(x) {
        x*dur_ICase + (1-x)*dur_IMild
      })

      rel_vacc <- out$odin_parameters$rel_infectiousness_vaccinated
      adjusted_eigens <- lapply(prop_susc, function(x) {

        unlist(lapply(seq_len(nrow(x)), function(y) {
          if(any(is.na(x[y,,]))) {
            return(NA)
          } else {
            Re(eigen(mixing_matrix*rowSums(x[y,,]*rel_vacc*relative_R0_by_age[[y]]))$values[1])
          }
        }))

      })

      ratios <-  (out$odin_parameters$beta_set * adjusted_eigens[[1]]) / out$parameters$R0

      # and patch NA gaps
      if(any(is.na(ratios))) {
        ratios[which(is.na(ratios))] <- ratios[which(!is.na(ratios))[1]]
      }

      ratios <- list(ratios)
    } else {
      prop_susc <- lapply(seq_len(dim(out$output)[3]), function(x) {
        t(t(out$output[seq_len(t_now), index$S, x])/pop)
      } )

      relative_R0_by_age <- prob_hosp*dur_ICase + (1-prob_hosp)*dur_IMild

      adjusted_eigens <- lapply(prop_susc, function(x) {

        unlist(lapply(seq_len(nrow(x)), function(y) {
          if(any(is.na(x[y,]))) {
            return(NA)
          } else {
            Re(eigen(mixing_matrix*x[y,]*relative_R0_by_age)$values[1])
          }
        }))

      })

      ratios <- list(
        out$odin_parameters$beta_set * adjusted_eigens[[1]]/out$parameters$R0
      )
    }
  } else {

    mixing_matrix <- squire:::process_contact_matrix_scaled_age(
      out$pmcmc_results$inputs$model_params$contact_matrix_set[[1]],
      out$pmcmc_results$inputs$model_params$population
    )

    dur_ICase <- out$parameters$dur_ICase
    dur_IMild <- out$parameters$dur_IMild
    prob_hosp <- out$parameters$prob_hosp

    # assertions
    squire:::assert_single_pos(dur_ICase, zero_allowed = FALSE)
    squire:::assert_single_pos(dur_IMild, zero_allowed = FALSE)
    squire:::assert_numeric(prob_hosp)
    squire:::assert_numeric(mixing_matrix)
    squire:::assert_square_matrix(mixing_matrix)
    squire:::assert_same_length(mixing_matrix[,1], prob_hosp)

    if(sum(is.na(prob_hosp)) > 0) {
      stop("prob_hosp must not contain NAs")
    }

    if(sum(is.na(mixing_matrix)) > 0) {
      stop("mixing_matrix must not contain NAs")
    }

    index <- squire:::odin_index(out$model)
    pop <- out$parameters$population

    if(is.null(max_date)) {
      max_date <- max(get_dates(out))
    }
    t_now <- which(as.Date(rownames(out$output)) == max_date)

    if(vaccine){
      # make the vaccine time changing args
      nrs <- t_now
      vei <- lapply(seq_len(nrow(out$pmcmc_results$inputs$model_params$vaccine_efficacy_infection)),
                    function(x) {out$pmcmc_results$inputs$model_params$vaccine_efficacy_infection[x,,]}
      )
      phl <- lapply(seq_len(nrow(out$pmcmc_results$inputs$model_params$prob_hosp)),
                    function(x) {out$pmcmc_results$inputs$model_params$prob_hosp[x,,]}
      )

      if(nrs > length(vei)) {
        vei_full <- c(rep(list(vei[[1]]),nrs - length(vei)), vei)
        phl_full <- c(rep(list(phl[[1]]),nrs - length(phl)), phl)
      } else {
        vei_full <- utils::tail(vei, nrs)
        phl_full <- utils::tail(phl, nrs)
      }

      # prop susceptible
      prop_susc <- lapply(seq_len(dim(out$output)[3]), function(x) {

        susceptible <- array(
          out$output[seq_len(t_now),index$S,x],
          dim=c(t_now, dim(index$S))
        )

        # We divide by the total population
        prop_susc <- sweep(susceptible, 2, pop, FUN='/')

        # We multiply by the effect of vaccines on onward infectiousness
        prop_susc <- vapply(
          seq_len(nrow(prop_susc)),
          FUN = function(i){prop_susc[i,,]*vei_full[[i]]},
          FUN.VALUE = prop_susc[1,,]
        )

        prop_susc <- aperm(prop_susc, c(3,1,2))

        return(prop_susc)
      } )

      relative_R0_by_age <- lapply(phl_full, function(x) {
        x*dur_ICase + (1-x)*dur_IMild
      })

      rel_vacc <- out$pmcmc_results$inputs$model_params$rel_infectiousness_vaccinated
      adjusted_eigens <- lapply(prop_susc, function(x) {

        unlist(lapply(seq_len(nrow(x)), function(y) {
          if(any(is.na(x[y,,]))) {
            return(NA)
          } else {
            Re(eigen(mixing_matrix*rowSums(x[y,,]*rel_vacc*relative_R0_by_age[[y]]))$values[1])
          }
        }))

      })

      betas <- lapply(out$replicate_parameters$R0, function(x) {
        squire:::beta_est(squire_model = out$pmcmc_results$inputs$squire_model,
                          model_params = out$pmcmc_results$inputs$model_params,
                          R0 = x)
      })

      ratios <- lapply(seq_along(betas), function(x) {
        (betas[[x]] * adjusted_eigens[[x]]) / out$replicate_parameters$R0[[x]]
      })

      # and patch NA gaps
      for (x in seq_along(ratios)) {
        if(any(is.na(ratios[[x]]))) {
          ratios[[x]][which(is.na(ratios[[x]]))] <- ratios[[x]][which(!is.na(ratios[[x]]))[1]]
        }
      }
    } else {
      prop_susc <- lapply(seq_len(dim(out$output)[3]), function(x) {
        t(t(out$output[seq_len(t_now), index$S, x])/pop)
      } )

      relative_R0_by_age <- prob_hosp*dur_ICase + (1-prob_hosp)*dur_IMild

      adjusted_eigens <- lapply(prop_susc, function(x) {

        unlist(lapply(seq_len(nrow(x)), function(y) {
          if(any(is.na(x[y,]))) {
            return(NA)
          } else {
            Re(eigen(mixing_matrix*x[y,]*relative_R0_by_age)$values[1])
          }
        }))

      })

      betas <- lapply(out$replicate_parameters$R0, function(x) {
        squire:::beta_est(squire_model = out$pmcmc_results$inputs$squire_model,
                          model_params = out$pmcmc_results$inputs$model_params,
                          R0 = x)
      })

      ratios <- lapply(seq_along(betas), function(x) {
        (betas[[x]] * adjusted_eigens[[x]]) / out$replicate_parameters$R0[[x]]
      })
    }
  }
  return(ratios)
}
#' An S3 method to calculate the proportion of the population immune
#' @inheritParams get_immunity_ratios
#'@export
get_immunity_ratios.rt_optimised <- function(out, max_date = NULL, vaccine = FALSE) {

  dates <- get_dates(out)

  #loop through each replicate
  purrr::map(seq_along(out$samples), function(sample_index){
    sample <- out$samples[[sample_index]]
    #remove initial infections since this doesn't matter (won't be regenerating simulations)
    sample$initial_infections <- NULL
    model_params <- setup_parameters(out$squire_model, append(out$parameters, sample))
    if(inherits(model_params$contact_matrix_set, "list")){
      mixing_matrix <- squire:::process_contact_matrix_scaled_age(
      model_params$contact_matrix_set[[1]],
      model_params$population
      )
    } else {
      mixing_matrix <- squire:::process_contact_matrix_scaled_age(
      model_params$contact_matrix_set,
      model_params$population
      )
    }

    dur_ICase <- 2/model_params$gamma_ICase
    dur_IMild <- 1/model_params$gamma_IMild

    # assertions
    squire:::assert_single_pos(dur_ICase, zero_allowed = FALSE)
    squire:::assert_single_pos(dur_IMild, zero_allowed = FALSE)
    squire:::assert_numeric(mixing_matrix)
    squire:::assert_square_matrix(mixing_matrix)

    if(!vaccine){
      prob_hosp <- model_params$prob_hosp_baseline
      squire:::assert_numeric(prob_hosp)
      if(sum(is.na(prob_hosp)) > 0) {
        stop("prob_hosp must not contain NAs")
      }
      squire:::assert_same_length(mixing_matrix[,1], prob_hosp)
    }

    if(sum(is.na(mixing_matrix)) > 0) {
      stop("mixing_matrix must not contain NAs")
    }

    index <- squire:::odin_index(out$model)
    pop <- model_params$population

    if(is.null(max_date)) {
      max_date <- max(dates)
    }

    t_now <- which(as.Date(rownames(out$output)) == max_date)
    if(vaccine & is.null(model_params$prob_hosp_baseline)){
      # make the vaccine time changing args
      nrs <- t_now
      vei <- purrr::map(seq_len(nrs), function(t){
        #find the change point equivalent
        eff_index <- max(which(model_params$tt_vaccine_efficacy_infection <= t))
        model_params$vaccine_efficacy_infection[eff_index,]
      })
      phl <- purrr::map(seq_len(nrs), function(t){
        #find the change point equivalent
        eff_index <- max(which(model_params$tt_vaccine_efficacy_disease <= t))
        outer(model_params$prob_hosp, model_params$vaccine_efficacy_disease[eff_index,], "*")
      })

      # prop susceptible
      susceptible <- array(
        out$output[seq_len(t_now),index$S,sample_index],
        dim=c(t_now, dim(index$S))
      )

      # We divide by the total population
      prop_susc <- sweep(susceptible, 2, pop, FUN='/')

      # We multiply by the effect of vaccines on onward infectiousness
      prop_susc <- vapply(
        seq_len(nrow(prop_susc)),
        FUN = function(i){sweep(prop_susc[i,,], 2, vei[[i]], FUN = "*")},
        FUN.VALUE = prop_susc[1,,]
      )
      prop_susc <- aperm(prop_susc, c(3,1,2))

      relative_R0_by_age <- lapply(phl, function(x) {
        x*dur_ICase + (1-x)*dur_IMild
      })

      rel_vacc <- model_params$rel_infectiousness_vaccinated

      adjusted_eigen <- unlist(lapply(seq_len(nrow(prop_susc)), function(y) {
        if(any(is.na(prop_susc[y,,]))) {
          return(NA)
        } else {
          Re(eigen(mixing_matrix*rowSums(prop_susc[y,,]* sweep(relative_R0_by_age[[y]], 2, rel_vacc, FUN = "*")))$values[1])
        }
      }))

      beta <- model_params$beta_set[1]

      ratio <- (beta * adjusted_eigen) / out$samples[[sample_index]]$R0[1]

      # and patch NA gaps
      if(any(is.na(ratio))) {
        ratio[which(is.na(ratio))] <- ratio[which(!is.na(ratio))[1]]
      }
    }
    else if(vaccine){
      # make the vaccine time changing args
      nrs <- t_now
      vei <- purrr::map(seq_len(nrs), function(t){
        #find the change point equivalent
        eff_index <- max(which(model_params$tt_vaccine_efficacy_infection <= t))
        model_params$vaccine_efficacy_infection[eff_index,,]
      })
      phl <- purrr::map(seq_len(nrs), function(t){
        #find the change point equivalent
        eff_index <- max(which(model_params$tt_vaccine_efficacy_disease <= t))
        model_params$prob_hosp[eff_index,,]
      })

      # prop susceptible
      susceptible <- array(
        out$output[seq_len(t_now),index$S,sample_index],
        dim=c(t_now, dim(index$S))
      )

      # We divide by the total population
      prop_susc <- sweep(susceptible, 2, pop, FUN='/')

      # We multiply by the effect of vaccines on onward infectiousness
      prop_susc <- vapply(
        seq_len(nrow(prop_susc)),
        FUN = function(i){prop_susc[i,,]*vei[[i]]},
        FUN.VALUE = prop_susc[1,,]
      )
      prop_susc <- aperm(prop_susc, c(3,1,2))

      relative_R0_by_age <- lapply(phl, function(x) {
        x*dur_ICase + (1-x)*dur_IMild
      })

      rel_vacc <- model_params$rel_infectiousness_vaccinated

      adjusted_eigen <- unlist(lapply(seq_len(nrow(prop_susc)), function(y) {
        if(any(is.na(prop_susc[y,,]))) {
          return(NA)
        } else {
          Re(eigen(mixing_matrix*rowSums(prop_susc[y,,]*rel_vacc*relative_R0_by_age[[y]]))$values[1])
        }
      }))

      beta <- model_params$beta_set[1]

      ratio <- (beta * adjusted_eigen) / out$samples[[sample_index]]$R0[1]

      # and patch NA gaps
      if(any(is.na(ratio))) {
        ratio[which(is.na(ratio))] <- ratio[which(!is.na(ratio))[1]]
      }
    } else {
      prop_susc <- t(t(out$output[seq_len(t_now), index$S, sample_index])/pop)

      relative_R0_by_age <- prob_hosp*dur_ICase + (1-prob_hosp)*dur_IMild

      adjusted_eigens <- unlist(lapply(seq_len(nrow(prop_susc)), function(y) {
        if(any(is.na(prop_susc[y,]))) {
          return(NA)
        } else {
          Re(eigen(mixing_matrix*prop_susc[y,]*relative_R0_by_age)$values[1])
        }
      }))

      beta <- model_params$beta_set[1]

      ratio <- (beta * adjusted_eigens) / out$samples[[sample_index]]$R0[1]
    }
    return(ratio)
  })
}
#' Plot the Rt trends for a country accounting for immunity
#' @param df model output
#' @param min_date Date to start the plot at
#' @param date_0 Current date
#' @param vjust unused
#' @param R0 plot R0
#' @param Rt plot Rt
#' @param Reff plot Reff
#' @export
country_immunity_plot <- function(df, min_date, date_0, vjust = -1.2, R0 = FALSE, Rt = FALSE, Reff = TRUE) {
  g1 <- ggplot2::ggplot(df %>% dplyr::filter(
    .data$date > min_date & .data$date <= as.Date(as.character(date_0 + as.POSIXlt(date_0)$wday + 1)))) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=12)) +
    ggplot2::xlab("") +
    ggplot2::ylab("Reff") +
    ggplot2::scale_x_date(breaks = "2 weeks",
                 limits = as.Date(c(as.character(min_date),
                                    as.character(date_0 + as.POSIXlt(date_0)$wday + 1))),
                 date_labels = "%d %b",
                 expand = c(0,0)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = "black")
    )

  if(Reff){
    g1 <- g1 +
      ggplot2::geom_ribbon(mapping = ggplot2::aes(x = .data$date, ymin=.data$Reff_min, ymax = .data$Reff_max, group = .data$iso3c), fill = "#96c4aa") +
      ggplot2::geom_ribbon(mapping = ggplot2::aes(x = .data$date, ymin = .data$Reff_q25, ymax = .data$Reff_q75, group = .data$iso3c), fill = "#48996b") +
      ggplot2::geom_line(mapping = ggplot2::aes(x = .data$date, y = .data$Reff_median), color = "#48996b")
  }

  if(R0) {
    g1 <- g1 +
      ggplot2::geom_ribbon(mapping = ggplot2::aes(x=.data$date, ymin=.data$R0_min, ymax = .data$R0_max, group = .data$iso3c), fill = "#8cbbca") +
      ggplot2::geom_ribbon(mapping = ggplot2::aes(x = .data$date, ymin = .data$R0_q25, ymax = .data$R0_q75, group = .data$iso3c), fill = "#3f8da7") +
      ggplot2::geom_hline(yintercept = df$R0_median[1], linetype = "dashed")
  }
  if(Rt) {
    g1 <- g1 +
      ggplot2::geom_ribbon(mapping = ggplot2::aes(x=.data$date, ymin=.data$Rt_min, ymax = .data$Rt_max, group = .data$iso3c), fill = "grey71", alpha = 0.25) +
      ggplot2::geom_ribbon(mapping = ggplot2::aes(x = .data$date, ymin = .data$Rt_q25, ymax = .data$Rt_q75, group = .data$iso3c), fill = "grey71", alpha = 0.55) +
      ggplot2::geom_line(mapping = ggplot2::aes(x = .data$date, y = .data$Rt_median), colour = "grey71", alpha = 0.55)
  }
  g1
}

#' Plot the Rt trends accounting for immunity
#' @param out model output
#' @param vaccine Should we account for vaccinations?
#' @param R0_plot plot R0
#' @param Rt_plot plot Rt
#'@export
rt_plot_immunity <- function(out, vaccine = TRUE, R0_plot = FALSE, Rt_plot = FALSE) {


  iso3c <- squire::get_population(out$parameters$country)$iso3c[1]

  date_0 <- get_data_end_date(out)

  # impact of immunity ratios
  ratios <- get_immunity_ratios(out, vaccine = vaccine)

  # get Rt values for each replicate
  rt <- get_Rt(out) %>% # add ratios
    dplyr::group_by(rep) %>%
    dplyr::mutate(
      ratios = ratios[[unique(.data$rep)]][seq_along(.data$Rt)],
      R0 = .data$Rt[1]*.data$ratios,
      Reff = .data$Rt*.data$ratios
    ) %>%
    dplyr::ungroup()

  rt$date <- as.Date(rt$date)

  new_rt_all <- rt %>%
    dplyr::group_by(rep) %>%
    dplyr::arrange(date)

  suppressMessages(sum_rt <- dplyr::group_by(new_rt_all, .data$date, .data$iso3c) %>%
                     dplyr::summarise(Rt_min = stats::quantile(.data$Rt, 0.025,na.rm=TRUE),
                               Rt_q25 = stats::quantile(.data$Rt, 0.25,na.rm=TRUE),
                               Rt_q75 = stats::quantile(.data$Rt, 0.75,na.rm=TRUE),
                               Rt_max = stats::quantile(.data$Rt, 0.975,na.rm=TRUE),
                               Rt_median = stats::median(.data$Rt,na.rm=TRUE),
                               Rt = mean(.data$Rt,na.rm=TRUE),
                               R0_min = stats::quantile(.data$R0, 0.025,na.rm=TRUE),
                               R0_q25 = stats::quantile(.data$R0, 0.25,na.rm=TRUE),
                               R0_q75 = stats::quantile(.data$R0, 0.75,na.rm=TRUE),
                               R0_max = stats::quantile(.data$R0, 0.975,na.rm=TRUE),
                               R0_median = stats::median(.data$R0),
                               R0 = mean(.data$R0),
                               Reff_min = stats::quantile(.data$Reff, 0.025,na.rm=TRUE),
                               Reff_q25 = stats::quantile(.data$Reff, 0.25,na.rm=TRUE),
                               Reff_q75 = stats::quantile(.data$Reff, 0.75,na.rm=TRUE),
                               Reff_max = stats::quantile(.data$Reff, 0.975,na.rm=TRUE),
                               Reff_median = stats::median(.data$Reff,na.rm=TRUE),
                               Reff = mean(.data$Reff,na.rm=TRUE)))

  min_date <- min(as.Date(out$inputs$start_date))

  res <- list("plot" = suppressWarnings(
    country_immunity_plot(sum_rt, min_date, date_0, R0 = R0_plot, Rt = Rt_plot, Reff = TRUE)
    ), "rts" = sum_rt)
  return(res)
}


#' Plot the seroprevalence from the model
#' Under some assumptions
#' @param res model output
#' @param sero_df Data to plot against
#'@export
sero_plot <- function(res, sero_df) {

  #compatibility with weekly fits
  date_0 <- get_data_end_date(res)

  # seroconversion data from brazeay report 34
  sero_det <- res$pmcmc_results$inputs$pars_obs$sero_det

  # additional_functions for rolling
  roll_func <- function(x, det) {
    l <- length(det)
    ret <- rep(0, length(x))
    for(i in seq_along(ret)) {
      to_sum <- utils::tail(x[seq_len(i)], length(det))
      ret[i] <- sum(rev(to_sum)*utils::head(det, length(to_sum)))
    }
    return(ret)
  }


  # get symptom onset data
  inf <- nimue_format(res, c("infections"), date_0 = date_0) %>%
    dplyr::rename(symptoms = "y") %>%
    dplyr::left_join(nimue_format(res, "S", date_0 = date_0),
              by = c("replicate", "t", "date")) %>%
    dplyr::rename(S = "y") %>%
    dplyr::select("replicate","t", "date", "S", "symptoms")

  inf <- inf %>% dplyr::group_by(.data$replicate) %>%
    dplyr::mutate(sero_positive = roll_func(.data$symptoms, .data$sero_det),
           sero_perc = .data$sero_positive/max(.data$S,na.rm = TRUE)) %>%
    dplyr::group_by(.data$date) %>%
    dplyr::summarise(sero_perc_med = stats::median(.data$sero_perc, na.rm=TRUE),
              sero_perc_min = stats::quantile(.data$sero_perc, 0.025, na.rm=TRUE),
              sero_perc_max = stats::quantile(.data$sero_perc, 0.975, na.rm=TRUE))

  gg <- ggplot2::ggplot(inf, ggplot2::aes(.data$date, .data$sero_perc_med, ymin = .data$sero_perc_min, ymax = .data$sero_perc_max)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(alpha = 0.2) +
    ggplot2::geom_point(ggplot2::aes(x = .data$date_start + (.data$date_end-.data$date_start)/2, y = .data$sero),
               sero_df, inherit.aes = FALSE) +
    ggplot2::geom_errorbar(ggplot2::aes(x = .data$date_start + (.data$date_end-.data$date_start)/2,
                      ymin = .data$sero_min, ymax = .data$sero_max),
                  sero_df, inherit.aes = FALSE, width = 0) +
    ggplot2::geom_errorbarh(ggplot2::aes(y = .data$sero, xmin = .data$date_start, xmax = .data$date_end),
                   sero_df, inherit.aes = FALSE, height = 0) +
    ggplot2::scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::theme_bw() + ggplot2::ylab("Seroprevalence") + ggplot2::xlab("")
  gg

  return(gg)

}

#' Plot the attack rate from the model
#' @param res model output
#'@export
ar_plot <- function(res) {

  #compatibility with weekly fits
  date_0 <- get_data_end_date(res)


  S_tot <- sum(get_parameters(res)$population)
  inf <- nimue_format(res, "infections", date_0 = date_0) %>%
    dplyr::mutate(infections = as.integer(.data$y)) %>%
    dplyr::select(replicate, "t", "date", "infections") %>%
    dplyr::group_by(replicate) %>%
    dplyr::mutate(infections = dplyr::lag(cumsum(tidyr::replace_na(.data$infections, 0)), 5, default = 0))


  g2 <- ggplot2::ggplot(inf, ggplot2::aes(.data$date, .data$infections/S_tot, group = .data$replicate)) + ggplot2::geom_line() +
    ggplot2::scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
    ggplot2::ylab("Attack Rate") + ggplot2::xlab("") + ggplot2::theme_bw()

  return(g2)
}

#' Plot the cumulative deaths from the model
#' @param res model output
#' @param extra_df data to plot if need
#'@export
cdp_plot <- function(res, extra_df = NULL) {

  date_0 <- get_data_end_date(res)

  #summarise deaths
  data <- get_data(res)
  data$date <- get_dates_greater(res) #assume reaches the cumulative sum at this
  #date works for both daily and weekly
  data$deaths <- cumsum(data$deaths)

  suppressWarnings(
    cdp <- plot(res, var_select = "D", date_0 = date_0, x_var = "date") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none", axis.title.x = ggplot2::element_blank()) +
      ggplot2::ylab("Cumulative Deaths") +
      ggplot2::scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
      ggplot2::xlab("") +
      ggplot2::geom_line(
        data = data,
        ggplot2::aes(
          x = .data$date, y = .data$deaths
        ),
        linetype = "dashed"
      )
  )

  if(!is.null(extra_df)){
    cdp <- cdp +
      ggplot2::geom_line(
        data = extra_df,
        ggplot2::aes(
          x = .data$date, y = cumsum(.data$deaths)
        ),
        alpha = 0.5,
        colour = "green"
      )+
      ggplot2::geom_ribbon(
        data = extra_df,
        ggplot2::aes(
          x = .data$date, ymin = cumsum(.data$bot), ymax = cumsum(.data$top)
        ),
        alpha = 0.25,
        fill = "green"
      )
  }

  cdp
}
#' Plot daily model deaths.
#'
#' Returns a ggplot2 object with median (95% CI) modelled deaths and the data it
#' has been fitted too.
#'
#' @param res A squire/nimue model object with generated fits (i.e. output of
#' generate_draws)
#' @return A ggplot2 object
#'@export
dp_plot <- function(res) {

  data <- get_data(res)

  data$date <- get_dates(res)

  suppressWarnings(
    dp <- plot(res, particle_fit = TRUE) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none", axis.title.x = ggplot2::element_blank()) +
      ggplot2::scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
      ggplot2::xlab("")
  )

  #compatibility with weekly fits
  if(any(c("rt_optimised", "excess_nimue_simulation") %in% class(res))){
    dp <- dp +
      ggplot2::geom_segment(data = get_data(res),
                   ggplot2::aes(x = .data$date_start, xend = .data$date_end,
                       y = .data$deaths/as.numeric(.data$date_end - .data$date_start),
                       yend = .data$deaths/as.numeric(.data$date_end - .data$date_start)),
                   linewidth = 1)
  } else {
    dp <- dp +
      ggplot2::geom_point(data = get_data(res),
                   ggplot2::aes(x = .data$date, y = .data$deaths),
                   size = 1)
  }

  dp

}
