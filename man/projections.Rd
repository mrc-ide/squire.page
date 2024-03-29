% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/projections.R
\name{projections}
\alias{projections}
\alias{projections.default}
\alias{projections.rt_optimised}
\title{Provide projections from calibrated simulations by changing RO, contact
matrices or bed availability.}
\usage{
projections(r, ...)

\method{projections}{default}(r, ...)

\method{projections}{rt_optimised}(
  r,
  time_period,
  R0 = NULL,
  R0_change = NULL,
  tt_R0 = NULL,
  contact_matrix_set = NULL,
  contact_matrix_set_change = NULL,
  tt_contact_matrix = NULL,
  hosp_bed_capacity = NULL,
  hosp_bed_capacity_change = NULL,
  tt_hosp_beds = NULL,
  ICU_bed_capacity = NULL,
  ICU_bed_capacity_change = NULL,
  tt_ICU_beds = NULL,
  to_be_run = TRUE,
  model_user_args = NULL,
  ...
)
}
\arguments{
\item{r}{Calibrated \code{{squire_simulation}} or \code{{rt_optimised}} object.}

\item{...}{any other parameters, only used to pass these parameters to the methods}

\item{time_period}{How many days is the projection. Default = NULL, which will
carry the projection forward from t = 0 in the calibration (i.e. the number
of days set in calibrate using forecast)}

\item{R0}{Numeric vector for R0 from t = 0 in the calibration.
E.g. \code{R0 = c(2, 1)}. Default = NULL, which will use \code{R0_change}
to alter R0 if provided.}

\item{R0_change}{Numeric vector for relative changes in R0 relative to the
final R0 used in the calibration (i.e. at t = 0 in the calibration)
E.g. \code{R0 = c(0.8, 0.5)}. Default = NULL, which will use \code{R0} to
parameterise changes in R0 if provided.}

\item{tt_R0}{Change time points for R0}

\item{contact_matrix_set}{Contact matrices used in simulation. Default =
NULL, which will use \code{contact_matrix_set_change} to alter the contact
matrix if provided.}

\item{contact_matrix_set_change}{Numeric vector for relative changes in the
contact matrix realtive to the final contact matrix used in the calibration
(i.e. at t = 0 in the calibration).
E.g. \code{contact_matrix_set_change = c(0.8, 0.5)}. Default = NULL, which
will use \code{contact_matrix_set} to parameterise changes in contact
matrices if if provided.}

\item{tt_contact_matrix}{Time change points for matrix change. Default = 0}

\item{hosp_bed_capacity}{Numeric vector for hospital bed capacity
from t = 0 in the calibration. Default = NULL, which will use
\code{hosp_bed_capacity_change} to alter hosp_bed_capacity if provided.}

\item{hosp_bed_capacity_change}{Numeric vector for relative changes in
hospital bed capacity relative to the final hospital bed capacity used in the
calibration (i.e. at t = 0 in the calibration).
E.g. \code{hosp_bed_capacity_change = c(0.8, 0.5)}. Default = NULL, which
will use \code{hosp_bed_capacity} to parameterise changes in hospital bed capacity
if provided.}

\item{tt_hosp_beds}{Change time points for hosp_bed_capacity}

\item{ICU_bed_capacity}{Numeric vector for ICU bed capacity
from t = 0 in the calibration. Default = NULL, which will use
\code{ICU_bed_capacity_change} to alter ICU_bed_capacity if provided.}

\item{ICU_bed_capacity_change}{Numeric vector for relative changes in
ICU bed capacity relative to the final ICU bed capacity used in the
calibration (i.e. at t = 0 in the calibration).
E.g. \code{ICU_bed_capacity_change = c(0.8, 0.5)}. Default = NULL, which
will use \code{ICU_bed_capacity} to parameterise changes in ICU bed capacity
if provided.}

\item{tt_ICU_beds}{Change time points for ICU_bed_capacity}

\item{to_be_run}{List of logicals for whether each replicate should be run.
Default = TRUE, which causes all replictes to be run.}

\item{model_user_args}{List of other parameters to be passed to the model to
be run. Default = NULL. An example would be:

\code{list(
  list(
  "prob_severe" = runif(17),
  "tt_dur_get_ox_survive" = c(0, 10),
  "gamma_get_ox_survive" = 0.2),
  list(
  "prob_severe" = runif(17),
  "tt_dur_get_mv_survive" = c(0, 5),
  "gamma_get_mv_survive" = 0.1)
  )}

The list should be the same length as the number of replicates in the
simulations. Each list element should then be a list with elements
named to match the arguments expected by the odin model with \code{r}.
Above would be suitable to set the model parameters for a simulation with
2 replicates. You do not have to have the same arguments in each list.}
}
\description{
This extends previous \code{projections} as you can pass in lists of each argument
that then get passed to each simulation replicate.
}
\details{
The user can specify changes to R0, contact matrices and bed
provision, which will come into effect from the current day in the calibration.
These changes can either set these to be specific values or change them
relative to their values in the original simulation. If no change is
requested, the simulation will use parameters chosen for the calibration run.
This extends previous versions of \code{projections} as you can now pass in
lists of each argument that then get passed to each simulation replicate.
}
