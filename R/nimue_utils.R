#' Nimue format that returns just compartments/summaries requested
#'
#' If a squire object is provided out it should be feed through to
#' \code{squire::format_output}. However, please note that not all variables can
#' be select for squire.
#'
#' @param out nimue_simulation or squire_simulation object
#' @param var_select Vector of compartment names, e.g. \code{c("S", "R")}. In
#'   addition a number of summary compartment can be requested. These include:
#' \itemize{
#'       \item{"deaths"}{ Daily Deaths }
#'       \item{"infections"}{ Daily Infections }
#'       \item{"hospital_occupancy"}{ Occupied Hospital Beds }
#'       \item{"ICU_occupancy"}{ Occupied ICU Beds }
#'       \item{"hospital_demand"}{ Required Hospital Beds }
#'       \item{"ICU_demand"}{ Required ICU Beds }
#'       \item{"hospital_incidence"}{ Incidence of hospitilisation }
#'       \item{"ICU_incidence"}{ Incidence of ICU admissions }
#'       \item{"hospitilisations"}{ Incidence of hospital+ICU admissions }
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
                "vaccinated_first_waned", "vaccinated_second_waned",
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
    if(any(c("vaccinated_first_dose", "vaccinated_second_dose", "vaccinated_first_waned", "vaccinated_second_waned",
             "first_doses_given", "second_doses_given", "booster_doses_given") %in% var_select)){
      warning(paste0(
        paste0(intersect(var_select, c("vaccinated_first_dose", "vaccinated_second_dose", "vaccinated_first_waned", "vaccinated_second_waned",
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

  # to match with squire uses
  if("hospital_incidence" %in% summs) {
    summs <- summs[-which(summs == "hospital_incidence")]
    hosp_inc_fix <- TRUE
  } else {
    hosp_inc_fix <- FALSE
  }

  # to match with squire uses
  if("ICU_incidence" %in% summs) {
    summs <- summs[-which(summs == "ICU_incidence")]
    ICU_inc_fix <- TRUE
  } else {
    ICU_inc_fix <- FALSE
  }

  # grab model outputs required
  if((length(comps) + length(summs)) != 0) {

    pd <- do.call(rbind, lapply(seq_len(dim(out$output)[3]), function(i) {
      format_squirepage(out, compartments = comps, summaries = summs, replicate = i, reduce_age = reduce_age)
    })) %>%
      dplyr::rename(y = .data$value)

    if(reduce_age){
      pd <- pd[,c("replicate", "compartment", "t", "y")]
    } else {
      pd <- pd[,c("replicate", "compartment", "age_group", "t", "y")]
    }

  } else {
    pd <- data.frame()
  }

  # add in hosp_inc
  if (hosp_inc_fix) {

    pd_hosp_inc <- do.call(rbind, lapply(seq_len(dim(out$output)[3]), function(i) {
      format_squirepage(out, compartments = "ICase2", summaries = character(0), replicate = i, reduce_age = FALSE)
    })) %>%
      dplyr::rename(y = .data$value) %>% dplyr::ungroup()

    prob_severe_age <- out$odin_parameters$prob_severe[as.numeric(pd_hosp_inc$age_group)]
    pd_hosp_inc$y <- out$odin_parameters$gamma_ICase * pd_hosp_inc$y * (1 - prob_severe_age)
    pd_hosp_inc$compartment <- "hospital_incidence"
    if(reduce_age) {
      pd_hosp_inc <- dplyr::group_by(pd_hosp_inc, .data$replicate, .data$compartment, .data$t) %>%
        dplyr::summarise(y = sum(.data$y))
    } else {
      pd_hosp_inc <- dplyr::group_by(pd_hosp_inc, .data$replicate, .data$compartment, .data$age_group,.data$t) %>%
        dplyr::summarise(y = sum(.data$y))
    }
    pd <- rbind(pd, pd_hosp_inc)
  }

  # add in ICU_inc
  if (ICU_inc_fix) {

    pd_ICU_inc <- do.call(rbind, lapply(seq_len(dim(out$output)[3]), function(i) {
      format_squirepage(out, compartments = "ICase2", summaries = character(0), replicate = i, reduce_age = FALSE)
    })) %>%
      dplyr::rename(y = .data$value) %>% dplyr::ungroup()

    prob_severe_age <- out$odin_parameters$prob_severe[as.numeric(pd_ICU_inc$age_group)]
    pd_ICU_inc$y <- out$odin_parameters$gamma_ICase * pd_ICU_inc$y * (prob_severe_age)
    pd_ICU_inc$compartment <- "ICU_incidence"
    if(reduce_age) {
      pd_ICU_inc <- dplyr::group_by(pd_ICU_inc, .data$replicate, .data$compartment, .data$t) %>%
        dplyr::summarise(y = sum(.data$y))
    } else {
      pd_ICU_inc <- dplyr::group_by(pd_ICU_inc, .data$replicate, .data$compartment, .data$age_group, .data$t) %>%
        dplyr::summarise(y = sum(.data$y))
    }
    pd <- rbind(pd, pd_ICU_inc)
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
