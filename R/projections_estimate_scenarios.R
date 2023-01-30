#' Calculates the Rt used in the short term scenarios
#' @param model_out a model object
#' @param forcast_days how far to calculate for
#'@export
get_future_Rt <- function(model_out, forcast_days){
  if("excess_nimue_simulation" %in% class(model_out)){
    stop("This function cannot be used with excess_nimue_simulation at the moment")
  }
  #cut overall Rt values up by Rt period and calculate all positive and negative
  #trends
  #get values used
  R0 <-  model_out$replicate_parameters[["R0"]]
  estimation_period <- model_out$pmcmc_results$inputs$Rt_args$Rt_rw_duration
  lockdown_date <- model_out$pmcmc_results$inputs$Rt_args$date_Meff_change
  #function to get the Rt from the effects and vice versa
  f_rt <- function(x, replicate){
    return(
      2*R0[replicate]*stats::plogis(x)
    )
  }
  inv_f_rt <- function(Rt, replicate){
    return(
      stats::qlogis(Rt/(2*R0[replicate]))
    )
  }
  #estimate Rt over time for each replicate in model fit
  rts <- get_Rt(model_out) %>%
    dplyr::group_by(.data$rep) %>%
    dplyr::mutate(#get the Rt in terms of untransformed effects
      Rt_effect = inv_f_rt(.data$Rt, unique(.data$rep)),
      #add an indicator if post lockdown effect
      post_lockdown = dplyr::if_else(
        .data$date > lockdown_date,
        TRUE,
        FALSE
      )
    )
  #now fit the model and estimate gradient
  fitting <- TRUE
  while(fitting){
    rts <- rts %>%
      dplyr::mutate(
        #split into fitting periods
        period = dplyr::if_else(
          .data$post_lockdown,
          ((.data$t - min(.data$t[.data$post_lockdown]))%/%estimation_period) + 1,
          0
        )
      ) %>%
      dplyr::group_by(.data$rep, .data$period) %>%
      #estimate gradient with lm
      dplyr::mutate(
        gradient =
          stats::lm(y~x, data = data.frame(x = .data$t,
                                           y = .data$Rt_effect))$coefficients[2]
      )
    #check we have at least on gradient of positive and negative
    pos_neg <- rts %>%
      dplyr::group_by(.data$rep, .data$period) %>%
      dplyr::filter(dplyr::n() > 2) %>%
      dplyr::group_by(.data$rep) %>%
      dplyr::filter(.data$period > 0) %>%
      dplyr::mutate(
        gradient_pos = .data$gradient >0
      ) %>%
      dplyr::select("rep", "gradient_pos") %>%
      unique() %>%
      dplyr::summarise(
        got_both = utils::head(.data$gradient_pos, 1) != utils::tail(.data$gradient_pos, 1)
      )
    if(!all(pos_neg$got_both)){
      #reduce estimation_period
      estimation_period <- estimation_period - 2
      #check if less than 4
      if(estimation_period < 4){
        warning("Unable to estimate scenario effects correctly for some replicates,
                removing the problematic replicates")
        #then we just drop the replicates that don't have both
        pos_neg_true <- pos_neg %>%
          dplyr::filter(.data$got_both)
        #check this is possible
        if(nrow(pos_neg_true) == 0){
          stop("No non-problematic trends avaiable")
        }
        rts <- rts %>%
          dplyr::filter(rep %in% pos_neg_true$rep)
        fitting <- FALSE
      }
    } else {
      fitting <- FALSE
    }
  }

  #calculate the optimistic/pessimistic gradients
  trends <- rts %>%
    dplyr::group_by(rep) %>%
    dplyr::filter(.data$period > 0) %>%
    dplyr::summarise(
      optimistic = stats::median(.data$gradient[.data$gradient < 0], 0.25),
      pessimistic = stats::median(.data$gradient[.data$gradient > 0], 0.75)
    )
  final_rts <- rts %>%
    dplyr::group_by(rep) %>%
    dplyr::filter(t == max(t))
  #for each replicate produce Rt values for the next forcast_days
  Rt_futures <- lapply(
    c("optimistic", "pessimistic"),
    function(trend){
      lapply(seq_len(nrow(trends)), function(row){
        f_rt(
          final_rts[row, ] %>% dplyr::pull(.data$Rt_effect) +
            trends[row, ] %>% dplyr::pull(trend)*seq(1, forcast_days, by = 1),
          replicate = row
        )
      })
    }
  )
  names(Rt_futures) <- c("optimistic", "pessimistic")
  return(Rt_futures)
}

#' Applies new Rt to the data
#' @param model_user_args list of values for each replicates
#' @param Rt Values to change Rt to
#' @param model_out model object
#'@export
update_Rt <- function(model_user_args, Rt, model_out){
  #get the Rt values for each replicate over time
  Rt_past <- get_Rt(model_out)
  #across each replicate/list in model_user_args
  for(i in seq_along(Rt)){
    #get previous beta and its times for the replicates
    beta_past <- squire:::beta_est(squire_model = model_out$pmcmc_results$inputs$squire_model,
                                   model_params = model_out$pmcmc_results$inputs$model_params,
                                   R0 = Rt_past %>% dplyr::filter(rep == i) %>% dplyr::pull(Rt)) %>%
      utils::tail(1)
    #calculate the beta equivalent to our target Rt
    beta <- squire:::beta_est(squire_model = model_out$pmcmc_results$inputs$squire_model,
                              model_params = model_out$pmcmc_results$inputs$model_params,
                              R0 = Rt[[i]])
    #attach to args
    model_user_args[[i]]$beta_set <- c(beta_past, beta)
    model_user_args[[i]]$tt_beta <- c(0, seq_along(beta))
  }
  return(model_user_args)
}
