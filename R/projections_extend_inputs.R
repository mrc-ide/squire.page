#'
#'@export
extend_vaccine_inputs <- function(vacc_inputs, time_period, out, end_date) {

  if("max_vaccine" %in% names(vacc_inputs)){
    # weekly mean vaccine distributions
    max_vaccine <- mean(utils::tail(vacc_inputs$max_vaccine,7))

    # assume at least 20% vaccinated by end of the year for meeting covax deadlines
    if(max_vaccine == 0) {
      end_of_year <- end_date
      lubridate::`month<-`(end_of_year, 12)
      lubridate::`day<-`(end_of_year, 31)
      max_vaccine <- round(sum(get_parameters(out)$pop)*0.2/as.integer(end_of_year-end_date))
    }
    tt_vaccine <- 0

    # efficacies best to just extend at the same rate
    vei <- vapply(seq_along(vacc_inputs$vaccine_efficacy_infection),
                  function(x) {
                    vacc_inputs$vaccine_efficacy_infection[[x]][1]
                  }, numeric(1))
    vei_new <- stats::predict(
      stats::lm(y~x, data.frame("x" = seq_along(vei), "y" = vei)),
      newdata = data.frame("x" = length(vei)+seq_len(time_period))
    )
    vei_new <- vapply(vei_new, min, numeric(1), 0.8)

    vaccine_efficacy_infection <- lapply(vei_new, rep, 17)
    tt_vaccine_efficacy_infection <- seq_along(vaccine_efficacy_infection)-1

    vaccine_efficacy_infection_odin_array <- nimue:::format_ve_i_for_odin(
      vaccine_efficacy_infection = vaccine_efficacy_infection,
      tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection
    )

    # efficacies best to just extend at the same rate
    ved <- vapply(seq_along(vacc_inputs$vaccine_efficacy_disease),
                  function(x) {
                    vacc_inputs$vaccine_efficacy_disease[[x]][1]
                  }, numeric(1))
    ved_new <- stats::predict(
      stats::lm(y~x, data.frame("x" = seq_along(ved), "y" = ved)),
      newdata = data.frame("x" = length(ved)+seq_len(time_period))
    )
    ved_new <- vapply(ved_new, min, numeric(1), 0.98)
    vaccine_efficacy_disease <- lapply(ved_new, rep, 17)
    tt_vaccine_efficacy_disease <- seq_along(vaccine_efficacy_disease)-1

    vaccine_efficacy_disease_odin_array <- nimue:::format_ve_d_for_odin(
      vaccine_efficacy_disease = vaccine_efficacy_disease,
      tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
      prob_hosp = out$parameters$prob_hosp
    )

    # combine into model_user_args for projections
    mua <- list(
      "max_vaccine" = max_vaccine,
      "tt_vaccine" = 0,
      "vaccine_efficacy_infection" = vaccine_efficacy_infection_odin_array,
      "tt_vaccine_efficacy_infection" = tt_vaccine_efficacy_infection,
      "prob_hosp" = vaccine_efficacy_disease_odin_array,
      "tt_vaccine_efficacy_disease" = tt_vaccine_efficacy_disease
    )
  } else if ("booster_doses" %in% names(vacc_inputs)){
    # weekly mean vaccine distributions
    first_doses <- mean(utils::tail(vacc_inputs$first_doses,7))
    second_doses <- mean(utils::tail(vacc_inputs$second_doses,7))
    booster_doses <- mean(utils::tail(vacc_inputs$booster_doses,7))

    #leave efficacies as they are

    # combine into model_user_args for projections
    mua <- list(
      "first_doses" = first_doses,
      "tt_first_doses" = 0,
      "second_doses" = second_doses,
      "tt_second_doses" = 0,
      "booster_doses" = booster_doses,
      "tt_booster_doses" = 0
    )
  } else {
    stop("Invalid vacc_inputs")
  }


  mua <- rep(list(mua), dim(out$output)[3])
  return(mua)

}
