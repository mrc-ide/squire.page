
#' Add long COVID outputs
#'
#' For use nimue_format
#'
#' @noRd
get_lc <- function(res,
                   mean_over_80_age = 82.5,
                   case_to_infection_ratio = 0.195,
                   reduce_age = FALSE) {

  # Long COVID probabilities per age group
  lc <- c(3.7,8.0,12.8,14.8,22.1,29.3,17.3)
  age <- c(7, 14.5, 21, 30, 42.5, 60, 78)
  lc_df <- data.frame(age = age, lc = lc)

  # Simple loess relationship to work out per 5-year age bin relationship
  squire_ages <- c(seq(2.5,77.5,5), mean_over_80_age)
  lc_squire <- suppressWarnings(
    stats::predict(stats::loess(lc~age, lc_df, span = 0.6),
            newdata =  data.frame("age"=squire_ages)
    ))
  lc_squire[1] <- lc_squire[2]-(lc_squire[3]-lc_squire[2])
  lc_squire[17] <- lc_squire[16]-(lc_squire[14]-lc_squire[15])

  # Format the
  infs <- nimue_format(res, "infections", reduce_age = FALSE)
  infs$cases <- infs$y*case_to_infection_ratio
  infs$lc <- infs$cases * lc_squire[as.numeric(infs$age_group)]/100
  infs$lc_prev <- zoo::rollsum(infs$lc, k = 12*7, fill = 0)
  infs$compartment <- "long_covid"

  if(reduce_age) {
    infs <- dplyr::group_by(infs, .data$replicate, .data$compartment, .data$t) %>%
      dplyr::summarise(y = sum(lc, na.rm=TRUE)) %>%
      dplyr::select(replicate, .data$compartment, .data$t, .data$y)
  } else {
    infs <- dplyr::group_by(infs, .data$replicate, .data$compartment, .data$age_group, .data$t) %>%
      dplyr::summarise(y = sum(lc, na.rm=TRUE)) %>%
      dplyr::select(replicate, .data$compartment, .data$age_group, .data$t, .data$y)
  }

  return(infs)
}

