#' Nimue format that returns just compartments/summaries requested
#'
#' If a squire object is provided out it should be feed through to
#' \code{squire::format_output}. However, please note that not all variables can
#' be selected for squire.
#'
#' To better define the summaries related to the healthcare pathway:
#' \describe{
#'   \item{demand}{The total number of requiring hospitalisation and ICU care
#'   at the given time, regardless of if they've received it}
#'   \item{occupancy}{The total number of those who received ICU or hospital care
#'   at the given time}
#'   \item{incidence}{The total number of new people who require hospital or ICU
#'   care}
#'   \item{hospitalisations}{Sum of new people recieving hospital or ICU care}
#' }
#' Note: this is confusing and seemingly inconsistent i.e. hospital_incidence
#' does not capture those recovering from ICU and hospitalisations is stated
#' to equal the sum of the incidences elsewhere.
#'
#' @param out nimue_simulation or squire_simulation object
#' @param var_select Vector of compartment names, e.g. \code{c("S", "R")}. In
#'   addition a number of summary compartment can be requested. These include:
#' \itemize{
#'       \item{"deaths"}{ Daily Deaths }
#'       \item{"infections"}{ Daily Infections }
#'       \item{"hospital_occupancy"}{ See description }
#'       \item{"ICU_occupancy"}{ See description }
#'       \item{"hospital_demand"}{ See description }
#'       \item{"ICU_demand"}{ See description }
#'       \item{"hospital_incidence"}{ See description }
#'       \item{"ICU_incidence"}{ See description }
#'       \item{"hospitalisations"}{ See description }
#'       \item{"vaccines", "unvaccinated", "vaccinated", "priorvaccinated"}{ Vaccine outputs }
#'       \item{"long_covid"}{ Long COVID estimates}
#'       }
#' @param reduce_age Collapse age-dimension, calculating the total in the
#'   compartment.
#' @param combine_compartments Collapse compartments of same type together
#'   (e.g. E1 and E2 -> E)
#' @param date_0 Date of time 0, if specified a date column will be added
#' @param mean_over_80_age Mean age of over 80 population for calculating long
#'   covid incidence. Default = 82.5
#' @param case_to_infection_ratio Ratio of infections to cases for long covid
#'   calculation. Default = 0.195 from UK ONS analyses
#' @return Formatted long data.frame
#' @examples
#' #generate some outputs for the afghanistan fit and pass through to format
#' #getting deaths and long covid estimates
#' \dontrun{
#' generate_draws(afg_fit) %>%
#'   nimue_format(var_select = c("deaths", "long_covid"))
#' }
#' @export
nimue_format <- function(out,
                         var_select = NULL,
                         reduce_age = TRUE,
                         combine_compartments = TRUE,
                         date_0 = NULL,
                         mean_over_80_age = 82.5,
                         case_to_infection_ratio = 0.195) {

  #if a particle fit we do a bit of adjusting so it plays nice with the squire functions
  if("rt_optimised" %in% class(out)){
    out$parameters$day_return <- TRUE
    out$parameters$replicates <- length(out$samples)
  }

  #a check for if a squire object then we use that instead
  if(!"nimue_simulation" %in% class(out)) {
    return(
      squire::format_output(out, var_select, reduce_age, combine_compartments, date_0)
    )
  }

  # work out what compartments are being plotted
  compartments = c("S", "E",
                   "IMild", "ICase", "IICU", "IHospital",
                   "IRec", "R", "D")
  summaries = c("N",
                "hospitalisations",
                "hospital_demand","hospital_occupancy",
                "ICU_demand", "ICU_occupancy",
                "vaccines", "unvaccinated", "vaccinated", "priorvaccinated",
                "vaccinated_first_dose", "vaccinated_second_dose",
                "vaccinated_second_waned", "vaccinated_booster_dose",
                "vaccinated_booster_waned",
                "first_doses_given", "second_doses_given", "booster_doses_given",
                "hospital_incidence", "ICU_incidence",
                "infections", "deaths")

  #check for correct types
  if("lmic_booster_nimue_simulation" %in% class(out)){
    if(any(c("vaccinated", "priorvaccinated", "vaccines") %in% var_select)){
      warning(paste0(paste0(
        intersect(c("vaccinated", "priorvaccinated", "vaccines"), var_select),
        collapse = ", "),
        " vaccinated and priorvaccinated cannot be output for this model type"))
    }
  } else {
    if(any(c("vaccinated_first_dose", "vaccinated_second_dose",
             "vaccinated_second_waned", "vaccinated_booster_dose",
             "vaccinated_booster_waned",
             "first_doses_given", "second_doses_given", "booster_doses_given") %in% var_select)){
      warning(paste0(
        paste0(intersect(var_select, c("vaccinated_first_dose", "vaccinated_second_dose",
                                       "vaccinated_second_waned", "vaccinated_booster_dose",
                                       "vaccinated_booster_waned",
                                       "first_doses_given", "second_doses_given", "booster_doses_given")), collapse = ", "),
        " cannot be output for this model type")
      )
    }
  }

  comps <- var_select[var_select %in% compartments]
  summs <- var_select[var_select %in% summaries]

  if(!combine_compartments){
    single_comps <- c("S", "IMild", "IICU", "D")
    comps_single <- intersect(comps, single_comps)
    comps_multi <- setdiff(comps, single_comps)
    comps <- c(
      comps_single,
      if(length(comps_multi) > 0){paste0(comps_multi, c(1, 2))}
      )
  }

  # to match with squire definition
  if("infections" %in% summs) {
    keep_E2 <- "E2" %in% comps
    if(!keep_E2){
      comps <- c(comps, "E2")
    }
    summs <- summs[-which(summs == "infections")]
    inf_fix <- TRUE
  } else {
    inf_fix <- FALSE
  }

  keep_ICase2 <- NULL #Catch for if we need to keep the compartment in tact
  # to match with squire uses
  if("hospital_incidence" %in% summs) {
    keep_ICase2 <- "ICase2" %in% comps
    if(!keep_ICase2){
      comps <- c(comps, "ICase2")
    }
    summs <- summs[-which(summs == "hospital_incidence")]
    hosp_inc_fix <- TRUE
  } else {
    hosp_inc_fix <- FALSE
  }

  # to match with squire uses
  if("ICU_incidence" %in% summs) {
    if(is.null(keep_ICase2)){
      keep_ICase2 <- "ICase2" %in% comps
      if(!keep_ICase2){
        comps <- c(comps, "ICase2")
      }
    }
    summs <- summs[-which(summs == "ICU_incidence")]
    ICU_inc_fix <- TRUE
  } else {
    ICU_inc_fix <- FALSE
  }

  # grab model outputs required
  if((length(comps) + length(summs)) != 0) {
    pd <- purrr::map_dfr(seq_len(dim(out$output)[3]), function(i) {
      format_squirepage(out, compartments = comps, summaries = summs, replicate = i, reduce_age = reduce_age & !any(c(ICU_inc_fix, hosp_inc_fix)))
    }) %>%
      dplyr::rename(y = .data$value)

    if(any(c(hosp_inc_fix, ICU_inc_fix))){
      #extract data for the fix
      pd_hosp_ICU_fix <- pd %>%
        dplyr::filter(.data$compartment == "ICase2")
      if(!keep_ICase2){
        pd <- pd %>%
          dplyr::filter(.data$compartment != "ICase2")
      }
      if(reduce_age){
        pd <- dplyr::group_by(pd, .data$replicate, .data$compartment, .data$t) %>%
          dplyr::summarise(y = sum(.data$y), .groups = "drop")
      }
    }

    if(reduce_age){
      pd <- pd[,c("replicate", "compartment", "t", "y")]
    } else {
      pd <- pd[,c("replicate", "compartment", "age_group", "t", "y")]
    }

  } else {
    pd <- data.frame()
  }

  #if we are doing any of these adjustments and its an Rt_optimise we need to
  #adjust for any changing parameters over time
  if(any(c(ICU_inc_fix, hosp_inc_fix, inf_fix)) & inherits(out, "rt_optimised")){
    get_severe <- any(c(ICU_inc_fix, hosp_inc_fix))
    rt_optimise_parameters <- purrr::map(seq_along(unique(pd$replicate)), function(rep){
      sample <- append(out$parameters, out$samples[[rep]])
      sample$initial_infections <- NULL
      sample$day_return <- NULL
      sample$replicates <- NULL
      pars <- setup_parameters(out$squire_model, sample)
      t <- pd %>% dplyr::filter(replicate == rep) %>% dplyr::pull(t) %>% unique()
      if(get_severe){
        prob_severe_age = t(t(matrix(pars$prob_severe, nrow = length(pars$prob_severe), ncol = length(t))) *
                              block_interpolate(t, pars$prob_severe_multiplier, pars$tt_prob_severe_multiplier))
        colnames(prob_severe_age) <- seq_len(ncol(prob_severe_age)) - 1
        prob_severe_age <- tibble::as_tibble(prob_severe_age, .name_repair = "minimal") %>%
          dplyr::mutate(age_group_num = seq_len(nrow(prob_severe_age))) %>%
          tidyr::pivot_longer(
            !.data$age_group_num,
            names_to = "t", values_to = "prob_severe"
          ) %>%
          dplyr::mutate(t = as.numeric(.data$t))
      } else {
        prob_severe_age <- NULL
      }
      list(
        gamma_E = pars$gamma_E,
        gamma_ICase = pars$gamma_ICase,
        prob_severe_age = prob_severe_age
      )
    })
    #split into components
    if(inf_fix){
      gamma_E_rt_optimise <- purrr::map_dfr(
        rt_optimise_parameters,
        ~c(gamma_E = .x$gamma_E),
        .id = "replicate") %>%
        dplyr::mutate(replicate = as.numeric(.data$replicate))
    }
    if(get_severe){
      prob_severe_rt_optimise <- purrr::map_dfr(rt_optimise_parameters, ~.x$prob_severe_age, .id = "replicate") %>%
        dplyr::mutate(replicate = as.numeric(.data$replicate))

      gamma_ICase_rt_optimise <- purrr::map_dfr(
        rt_optimise_parameters,
        ~c(gamma_ICase = .x$gamma_ICase),
        .id = "replicate") %>%
        dplyr::mutate(replicate = as.numeric(.data$replicate))
    }
    rm(get_severe, rt_optimise_parameters)
  }

  # fix the infection
  if (inf_fix) {
    if(inherits(out, "rt_optimised")){
      new_vals <- pd %>%
        dplyr::filter(.data$compartment == "E2") %>%
        dplyr::left_join(
          gamma_E_rt_optimise,
          by = "replicate"
        ) %>%
        dplyr::mutate(
          y = .data$y *.data$gamma_E
        ) %>%
        dplyr::pull(.data$y)
      rm(gamma_E_rt_optimise)
    } else {
      new_vals <- pd$y[pd$compartment == "E2"] * out$odin_parameters$gamma_E
    }
    if(keep_E2){
      pd <- pd %>% rbind(
        pd %>%
          dplyr::filter(
            .data$compartment == "E2"
          ) %>%
          dplyr::mutate(
            y = new_vals,
            comparment = "infections"
          )
      )
    } else {
      pd$y[pd$compartment == "E2"] <- new_vals
      pd$compartment <- as.character(pd$compartment)
      pd$compartment[pd$compartment == "E2"] <- "infections"
      pd$compartment <- as.factor(pd$compartment)
    }
    rm(new_vals)
  }

  # add in hosp_inc
  if (hosp_inc_fix) {
    pd_hosp_ICU_fix
    if(inherits(out, "rt_optimised")){
      df <- pd_hosp_ICU_fix %>%
        dplyr::left_join(
          gamma_ICase_rt_optimise, by = "replicate"
        ) %>%
        dplyr::mutate(age_group_num = as.numeric(.data$age_group)) %>%
        dplyr::left_join(
          prob_severe_rt_optimise, by = c("replicate", "age_group_num", "t")
        ) %>%
        dplyr::mutate(
          y = .data$gamma_ICase * .data$y * (1 - .data$prob_severe),
          compartment = "hospital_incidence"
        )
      if(!ICU_inc_fix){
        rm(gamma_ICase_rt_optimise, prob_severe_rt_optimise)
      }
    } else {
      prob_severe_age <- out$odin_parameters$prob_severe[as.numeric(pd_hosp_ICU_fix$age_group)]
      df <- pd_hosp_ICU_fix
      df$y <- out$odin_parameters$gamma_ICase * df$y * (1 - (prob_severe_age *
        #add adjustment for prob severe multiplier (for Rt optimise this is already in prob_severe)
        block_interpolate(df$t, out$odin_parameters$prob_severe_multiplier, out$odin_parameters$tt_prob_severe_multiplier)
      ))
      df$compartment <- "hospital_incidence"
    }
    if(!ICU_inc_fix){
      rm(pd_hosp_ICU_fix)
    }
    if(reduce_age) {
      df <- dplyr::group_by(df, .data$replicate, .data$compartment, .data$t) %>%
        dplyr::summarise(y = sum(.data$y), .groups = "drop")
    } else {
      df <- dplyr::group_by(df, .data$replicate, .data$compartment, .data$age_group,.data$t) %>%
        dplyr::summarise(y = sum(.data$y), .groups = "drop")
    }
    pd <- rbind(pd, df)
  }

  # add in ICU_inc
  if (ICU_inc_fix) {
    pd_hosp_ICU_fix
    if(inherits(out, "rt_optimised")){
      df <- pd_hosp_ICU_fix %>%
        dplyr::left_join(
          gamma_ICase_rt_optimise, by = "replicate"
        ) %>%
        dplyr::mutate(age_group_num = as.numeric(.data$age_group)) %>%
        dplyr::left_join(
          prob_severe_rt_optimise, by = c("replicate", "age_group_num", "t")
        ) %>%
        dplyr::mutate(
          y = .data$gamma_ICase * .data$y * .data$prob_severe,
          compartment = "ICU_incidence"
        )
      rm(gamma_ICase_rt_optimise, prob_severe_rt_optimise)
    } else {
      prob_severe_age <- out$odin_parameters$prob_severe[as.numeric(pd_hosp_ICU_fix$age_group)]
      df <- pd_hosp_ICU_fix
      df$y <- out$odin_parameters$gamma_ICase * df$y * prob_severe_age *
        #add adjustment for prob severe multiplier (for Rt optimise this is already in prob_severe)
        block_interpolate(df$t, out$odin_parameters$prob_severe_multiplier, out$odin_parameters$tt_prob_severe_multiplier)
      df$compartment <- "ICU_incidence"
    }
    rm(pd_hosp_ICU_fix)
    if(reduce_age) {
      df <- dplyr::group_by(df, .data$replicate, .data$compartment, .data$t) %>%
        dplyr::summarise(y = sum(.data$y), .groups = "drop")
    } else {
      df <- dplyr::group_by(df, .data$replicate, .data$compartment, .data$age_group,.data$t) %>%
        dplyr::summarise(y = sum(.data$y), .groups = "drop")
    }
    pd <- rbind(pd, df)
  }

  # add in long covid too
   if("long_covid" %in% var_select) {
     lc <- get_lc(out, mean_over_80_age, case_to_infection_ratio, reduce_age)
     pd <- rbind(pd, lc)
   }

  # replacing time with date if date_0 is provided
  if(!is.null(date_0)){
    pd$date <- as.Date(pd$t + as.Date(date_0),
                       format = "%Y-%m-%d")
  }
  return(pd)
}

#' Format vaccine model output
#'
#' Take raw odin vaccine model output and formats in long format with the option to select
#' variables and summarise over age groups. Output variables are ordered as in argument ordering.
#'
#' @param x squire_simulation object
#' @param compartments Vector of compartment names, e.g. \code{c("S", "R")}, or sub-compartment names, e.g. \code{c("S", "E1", "E2")}
#' @param summaries Vector of summary names, which may be:
#' \itemize{
#'       \item{"deaths"}{ Deaths per day }
#'       \item{"infections"}{ Infections per day. New infections (note this is currently a slightly different definitionto the main Squire mode)}
#'       \item{"hospitilisations"}{ Hospitalisations per day (Note this takes into account hospital capacity)}
#'       \item{"hospital_occupancy"}{ Occupied Hospital Beds }
#'       \item{"ICU_occupancy"}{ Occupied ICU Beds }
#'       \item{"hospital_demand}{ Required Hospital Beds }
#'       \item{"ICU_demand}{ Required ICU Beds }
#'       \item{"vaccinated"}{ Vaccines administered per day}
#'       }
#' @param reduce_age Collapse age-dimension, calculating the total in the
#'   compartment.
#' @param date_0 Date of time 0 (e.g. "2020-03-01"), if specified a date column will be added
#' @param replicate Which replicate is being formatted. Default = 1
#'
#' @return Formatted long data.frame
#' @noRd
format_squirepage <- function(x,
                   compartments,
                   summaries,
                   reduce_age = TRUE,
                   date_0 = NULL,
                   replicate = 1){

  # Arg checks
  squire:::assert_custom_class(x, "nimue_simulation")
  squire:::assert_logical(reduce_age)

  # Standardise output dimensions
  if(length(dim(x$output)) == 4){
    x$output <- abind::adrop(x$output, drop = c(FALSE, FALSE, FALSE, TRUE))
  }

  # Get columns indices of variables
  index <- squire:::odin_index(x$model)
  if(!all(compartments %in% names(index))){
    stop("Some compartments specified not output by model")
  }

  # Extract time
  time <- x$output[,index$time, replicate]

  output <- nimue:::format_internal(x = x, compartments = compartments, summaries = summaries,
                            reduce_age = reduce_age, index = index,
                            time = time, replicate = replicate)

  # Set levels (order) of output variables
  output$compartment <- factor(output$compartment, levels = c(compartments, summaries))

  # Add date
  if(!is.null(date_0)){
    squire:::assert_date(date_0)
    output$date <- as.Date(output$t + as.Date(date_0),
                           format = "%Y-%m-%d")
  }

  # Add age-groups if present
  if("age_index" %in% names(output)){
    ag <- c(paste0(seq(0, 75, 5), "-", seq(5, 80, 5)), "80+")
    output$age_group = factor(ag[output$age_index], levels = ag)
    output <- output  %>%
      dplyr::select(-.data$age_index)
  }

  return(output)
}
