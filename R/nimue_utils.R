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
  if(!"tt_vaccine" %in% out$model$.__enclos_env__$private$user) {
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
                "hospital_incidence", "ICU_incidence",
                "infections", "deaths")

  comps <- var_select[var_select %in% compartments]
  summs <- var_select[var_select %in% summaries]

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
      nimue::format(out, compartments = comps, summaries = summs, replicate = i, reduce_age = reduce_age)
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
      nimue::format(out, compartments = "ICase2", summaries = character(0), replicate = i, reduce_age = FALSE)
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
      nimue::format(out, compartments = "ICase2", summaries = character(0), replicate = i, reduce_age = FALSE)
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