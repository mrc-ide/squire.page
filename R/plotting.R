#'
#'@export
compare_adjustment_plot <- function(out, grad_dur, extend_past = 10){
  #get the model infections
  plotting_data <- nimue_format(out, var_select = "infections") %>%
    dplyr::filter(t >= (max(.data$t, na.rm = TRUE) - grad_dur - extend_past)) %>%
    dplyr::rename(model_infections = .data$y) %>%
    dplyr::select(replicate, .data$t, .data$model_infections) %>%
    dplyr::mutate(
      reporting = .data$t < max(.data$t, na.rm = TRUE) - grad_dur
    ) %>%
    #add the real data
    dplyr::left_join(
      out$pmcmc_results$inputs$data %>%
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
    ggplot2::geom_vline(ggplot2::aes(xintercept = max(.data$date)-grad_dur),
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

#'
#'@export
get_immunity_ratios <- function(out, max_date = NULL) {

  out$pmcmc_results$inputs$data$date <- get_dates(out)

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
    max_date <- max(out$pmcmc_results$inputs$data$date)
  }
  t_now <- which(as.Date(rownames(out$output)) == max_date)
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

  return(ratios)
}
#'
#'@export
get_immunity_ratios_vaccine <- function(out, max_date = NULL) {

  out$pmcmc_results$inputs$data$date <- get_dates(out)

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
    max_date <- max(out$pmcmc_results$inputs$data$date)
  }
  t_now <- which(as.Date(rownames(out$output)) == max_date)

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

  return(ratios)
}
#'
#'@export
country_immunity_plot <- function(df, min_date, date_0, vjust = -1.2, R0 = FALSE, Rt = FALSE, Reff = TRUE) {
  g1 <- ggplot2::ggplot(df %>% dplyr::filter(
    .data$date > min_date & .data$date <= as.Date(as.character(date_0+as.numeric(lubridate::wday(date_0)))))) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=12)) +
    ggplot2::xlab("") +
    ggplot2::ylab("Reff") +
    ggplot2::scale_x_date(breaks = "2 weeks",
                 limits = as.Date(c(as.character(min_date),
                                    as.character(date_0+as.numeric(lubridate::wday(date_0))))),
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

#'
#'@export
rt_plot_immunity <- function(out, vaccine = TRUE, R0_plot = FALSE, Rt_plot = FALSE) {


  iso3c <- squire::get_population(out$parameters$country)$iso3c[1]

  date_0 <- get_data_end_date(out)

  # impact of immunity ratios
  if(vaccine){
    ratios <- get_immunity_ratios_vaccine(out)
  } else {
    ratios <- get_immunity_ratios(out)
  }

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

  min_date <- min(as.Date(out$replicate_parameters$start_date))

  res <- list("plot" = suppressWarnings(
    country_immunity_plot(sum_rt, min_date, date_0, R0 = R0_plot, Rt = Rt_plot, Reff = TRUE)
    ), "rts" = sum_rt)
  return(res)
}

#'
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
    dplyr::rename(symptoms = .data$y) %>%
    dplyr::left_join(nimue_format(res, "S", date_0 = date_0),
              by = c("replicate", "t", "date")) %>%
    dplyr::rename(S = .data$y) %>%
    dplyr::select(.data$replicate,.data$t, .data$date, .data$S, .data$symptoms)

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
#'
#'@export
ar_plot <- function(res) {

  #compatibility with weekly fits
  date_0 <- get_data_end_date(res)


  S_tot <- sum(res$pmcmc_results$inputs$model_params$population)
  inf <- nimue_format(res, "infections", date_0 = date_0) %>%
    dplyr::mutate(infections = as.integer(.data$y)) %>%
    dplyr::select(replicate, .data$t, .data$date, .data$infections) %>%
    dplyr::group_by(replicate) %>%
    dplyr::mutate(infections = dplyr::lag(cumsum(tidyr::replace_na(.data$infections, 0)), 5, default = 0))


  g2 <- ggplot2::ggplot(inf, ggplot2::aes(.data$date, .data$infections/S_tot, group = .data$replicate)) + ggplot2::geom_line() +
    ggplot2::scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
    ggplot2::ylab("Attack Rate") + ggplot2::xlab("") + ggplot2::theme_bw()

  return(g2)
}

#'
#'@export
cdp_plot <- function(res, extra_df = NULL) {

  date_0 <- get_data_end_date(res)

  #summarise deaths
  data <- res$pmcmc_results$inputs$data
  data$date <- get_dates_greater(res) #assume reaches the cumulative sum at this
  #date works for both daily and weekly
  data$deaths <- cumsum(data$deaths)

  suppressWarnings(
    cdp <- plot(res, "D", date_0 = date_0, x_var = "date") +
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
#'
#'@export
dp_plot <- function(res) {

  res$pmcmc_results$inputs$data$date <- get_dates(res)

  suppressWarnings(
    dp <- plot(res, particle_fit = TRUE) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none", axis.title.x = ggplot2::element_blank()) +
      ggplot2::scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
      ggplot2::xlab("")
  )

  #compatibility with weekly fits
  if("week_start" %in% names(res$pmcmc_results$inputs$data)){
    dp <- dp +
      ggplot2::geom_segment(data = res$pmcmc_results$inputs$data,
                   ggplot2::aes(x = .data$week_start, xend = .data$week_end,
                       y = .data$deaths/7, yend = .data$deaths/7),
                   size = 1)
    dp$layers[[5]] <- NULL
  } else {
    dp <- dp +
      ggplot2::geom_point(data = res$pmcmc_results$inputs$data,
                   ggplot2::aes(x = .data$date, y = .data$deaths),
                   size = 1)
  }

  dp

}
