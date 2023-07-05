################################################################################
### Vaccine model: deterministic ###############################################
################################################################################

################################################################################

## Time output to line up with squire fitting infrastructure
time <- t
output(time) <- TRUE


### S: susceptible #############################################################
dim(S) <- c(17, 8)

S_0[, ] <- user()
dim(S_0) <- c(17, 8)
initial(S[, ]) <- S_0[i, j]

#ISSUE with these empty rates might actually be quite slow see if it speeds up when changed
deriv(S[, 1]) <- (gamma_R_t * R2[i, j]) - (lambda[i] * vaccine_efficacy_infection_t[i, j] * S[i, j]) - gamma_vaccine_t[j] * S[i, j] + vaccinations_S[i,j]
deriv(S[, 2:8]) <- (gamma_R_t * R2[i, j]) - (lambda[i] * vaccine_efficacy_infection_t[i, j] * S[i, j]) - gamma_vaccine_t[j] * S[i, j] + vaccinations_S[i,j] + gamma_vaccine_t[j - 1] * S[i, j - 1]

vaccinations_S[, 1] <-  - primary_first[i] * S[i, j]
vaccinations_S[, 2] <-  primary_first[i] * S[i, j - 1] - primary_second[i] * S[i, j]
vaccinations_S[, 3] <-  primary_second[i] * S[i, j - 1] - booster_first[i] * S[i, j]
vaccinations_S[, 4] <-  - booster_first[i] * S[i, j]
vaccinations_S[, 5] <-  - booster_first[i] * S[i, j]
vaccinations_S[, 6] <-  booster_first[i] * sum(S[i, 3:5]) + booster_second[i] * sum(S[i, (j + 1):8])
vaccinations_S[, 7:8] <- -booster_second[i] * S[i, j]
dim(vaccinations_S) <- c(17, 8)

################################################################################

### E (E1 & E2): Latent ########################################################
dim(E1) <- c(17, 8)
dim(E2) <- c(17, 8)

E1_0[, ] <- user()
dim(E1_0) <- c(17, 8)
initial(E1[, ]) <- E1_0[i, j]

E2_0[, ] <- user()
dim(E2_0) <- c(17, 8)
initial(E2[, ]) <- E2_0[i, j]

gamma_E[] <- user() # rate of progression through latent infection
dim(gamma_E) <- user()
tt_dur_E[] <- user()
dim(tt_dur_E) <- length(gamma_E)
gamma_E_t <- interpolate(tt_dur_E, gamma_E, "constant")

deriv(E1[, 1]) <- (lambda[i] * vaccine_efficacy_infection_t[i, j] * S[i, j]) - (gamma_E_t * E1[i, j]) - gamma_vaccine_t[j] * E1[i, j] + vaccinations_E1[i,j]
deriv(E1[, 2:8]) <- (lambda[i] * vaccine_efficacy_infection_t[i, j] * S[i, j]) - (gamma_E_t * E1[i, j]) - gamma_vaccine_t[j] * E1[i, j] + vaccinations_E1[i,j] + gamma_vaccine_t[j - 1] * E1[i, j - 1]

vaccinations_E1[, 1] <-  - primary_first[i] * E1[i, j]
vaccinations_E1[, 2] <-  primary_first[i] * E1[i, j - 1] - primary_second[i] * E1[i, j]
vaccinations_E1[, 3] <-  primary_second[i] * E1[i, j - 1] - booster_first[i] * E1[i, j]
vaccinations_E1[, 4] <-  - booster_first[i] * E1[i, j]
vaccinations_E1[, 5] <-  - booster_first[i] * E1[i, j]
vaccinations_E1[, 6] <-  booster_first[i] * sum(E1[i, 3:5]) + booster_second[i] * sum(E1[i, (j + 1):8])
vaccinations_E1[, 7:8] <- -booster_second[i] * E1[i, j]
dim(vaccinations_E1) <- c(17, 8)

deriv(E2[, 1]) <- (gamma_E_t * E1[i, j]) - (gamma_E_t * E2[i, j]) - gamma_vaccine_t[j] * E2[i, j] + vaccinations_E2[i,j]
deriv(E2[, 2:8]) <- (gamma_E_t * E1[i, j]) - (gamma_E_t * E2[i, j]) - gamma_vaccine_t[j] * E2[i, j] + vaccinations_E2[i,j] + gamma_vaccine_t[j - 1] * E2[i, j - 1]

vaccinations_E2[, 1] <-  - primary_first[i] * E2[i, j]
vaccinations_E2[, 2] <-  primary_first[i] * E2[i, j - 1] - primary_second[i] * E2[i, j]
vaccinations_E2[, 3] <-  primary_second[i] * E2[i, j - 1] - booster_first[i] * E2[i, j]
vaccinations_E2[, 4] <-  - booster_first[i] * E2[i, j]
vaccinations_E2[, 5] <-  - booster_first[i] * E2[i, j]
vaccinations_E2[, 6] <-  booster_first[i] * sum(E2[i, 3:5]) + booster_second[i] * sum(E2[i, (j + 1):8])
vaccinations_E2[, 7:8] <- -booster_second[i] * E2[i, j]
dim(vaccinations_E2) <- c(17, 8)

output(E[]) <- sum(E1[i, ]) + sum(E2[i, ])
dim(E) <- 17
################################################################################

### IMild: Unhospitalised infection ############################################
dim(IMild) <- c(17, 8)

IMild_0[, ] <- user()
dim(IMild_0) <- c(17, 8)
initial(IMild[, ]) <- IMild_0[i, j]

gamma_IMild[] <- user() # rate of progression from mild infection to recovery
dim(gamma_IMild) <- user()
tt_dur_IMild[] <- user()
dim(tt_dur_IMild) <- length(gamma_IMild)
gamma_IMild_t <- interpolate(tt_dur_IMild, gamma_IMild, "constant")

deriv(IMild[, 1]) <- (gamma_E_t * E2[i, j] * (1 - prob_hosp_t_mult[i, j])) - (gamma_IMild_t * IMild[i, j]) - gamma_vaccine_t[j] * IMild[i, j]
deriv(IMild[, 2:8]) <- (gamma_E_t * E2[i, j] * (1 - prob_hosp_t_mult[i, j])) - (gamma_IMild_t * IMild[i, j]) - gamma_vaccine_t[j] * IMild[i, j] + gamma_vaccine_t[j - 1] * IMild[i, j - 1]

output(IMildout[]) <- sum(IMild[i, ])
dim(IMildout) <- 17
################################################################################

### R: (R1 & R2): Recovered ####################################################
dim(R1) <- c(17, 8)
dim(R2) <- c(17, 8)

R1_0[, ] <- user()
dim(R1_0) <- c(17, 8)
initial(R1[, ]) <- R1_0[i, j]

R2_0[, ] <- user()
dim(R2_0) <- c(17, 8)
initial(R2[, ]) <- R2_0[i, j]

# Interpolation for dur_R
gamma_R_t <- interpolate(tt_dur_R, gamma_R, "constant")
tt_dur_R[] <- user()
dim(tt_dur_R) <- user()
dim(gamma_R) <- length(tt_dur_R)
gamma_R[] <- user() # rate of progression through recovered compartment (loss of naturally acquired immunity)

deriv(R1[, 1]) <- (gamma_rec * IRec2[i, j]) + (gamma_IMild_t * IMild[i, j]) + (gamma_get_ox_survive_t * IOxGetLive2[i, j]) + (gamma_not_get_ox_survive * IOxNotGetLive2[i, j]) + (gamma_not_get_mv_survive * IMVNotGetLive2[i, j]) - (gamma_R_t * R1[i, j]) - gamma_vaccine_t[j] * R1[i, j] + vaccinations_R1[i,j]
deriv(R1[, 2:8]) <- (gamma_rec * IRec2[i, j]) + (gamma_IMild_t * IMild[i, j]) + (gamma_get_ox_survive_t * IOxGetLive2[i, j]) + (gamma_not_get_ox_survive * IOxNotGetLive2[i, j]) + (gamma_not_get_mv_survive * IMVNotGetLive2[i, j]) - (gamma_R_t * R1[i, j]) - gamma_vaccine_t[j] * R1[i, j] + vaccinations_R1[i,j] + gamma_vaccine_t[j - 1] * R1[i, j - 1]

vaccinations_R1[, 1] <-  - primary_first[i] * R1[i, j]
vaccinations_R1[, 2] <-  primary_first[i] * R1[i, j - 1] - primary_second[i] * R1[i, j]
vaccinations_R1[, 3] <-  primary_second[i] * R1[i, j - 1] - booster_first[i] * R1[i, j]
vaccinations_R1[, 4] <-  - booster_first[i] * R1[i, j]
vaccinations_R1[, 5] <-  - booster_first[i] * R1[i, j]
vaccinations_R1[, 6] <-  booster_first[i] * sum(R1[i, 3:5]) + booster_second[i] * sum(R1[i, (j + 1):8])
vaccinations_R1[, 7:8] <- -booster_second[i] * R1[i, j]
dim(vaccinations_R1) <- c(17, 8)

deriv(R2[, 1]) <- (gamma_R_t * R1[i, j]) - (gamma_R_t * R2[i, j]) - gamma_vaccine_t[j] * R2[i, j] + vaccinations_R2[i,j]
deriv(R2[, 2:8]) <- (gamma_R_t * R1[i, j]) - (gamma_R_t * R2[i, j]) - gamma_vaccine_t[j] * R2[i, j] + vaccinations_R2[i,j] + gamma_vaccine_t[j - 1] * R2[i, j - 1]

vaccinations_R2[, 1] <-  - primary_first[i] * R2[i, j]
vaccinations_R2[, 2] <-  primary_first[i] * R2[i, j - 1] - primary_second[i] * R2[i, j]
vaccinations_R2[, 3] <-  primary_second[i] * R2[i, j - 1] - booster_first[i] * R2[i, j]
vaccinations_R2[, 4] <-  - booster_first[i] * R2[i, j]
vaccinations_R2[, 5] <-  - booster_first[i] * R2[i, j]
vaccinations_R2[, 6] <-  booster_first[i] * sum(R2[i, 3:5]) + booster_second[i] * sum(R2[i, (j + 1):8])
vaccinations_R2[, 7:8] <- -booster_second[i] * R2[i, j]
dim(vaccinations_R2) <- c(17, 8)

output(R[]) <- sum(R1[i, ]) + sum(R2[i, ])
dim(R) <- 17
################################################################################

### ICase (ICase1 & ICase2): To-be hospitalised infection ######################
dim(ICase1) <- c(17, 8)
dim(ICase2) <- c(17, 8)

ICase1_0[, ] <- user()
dim(ICase1_0) <- c(17, 8)
initial(ICase1[, ]) <- ICase1_0[i, j]

ICase2_0[, ] <- user()
dim(ICase2_0) <- c(17, 8)
initial(ICase2[, ]) <- ICase2_0[i, j]

gamma_ICase[] <- user() # rate of progression from symptom onset to requiring hospitalisation
dim(gamma_ICase) <- user()
tt_dur_ICase[] <- user()
dim(tt_dur_ICase) <- length(gamma_ICase)
gamma_ICase_t <- interpolate(tt_dur_ICase, gamma_ICase, "constant")

deriv(ICase1[, 1]) <- (gamma_E_t * E2[i, j] * prob_hosp_t_mult[i, j]) - (gamma_ICase_t * ICase1[i, j]) - gamma_vaccine_t[j] * ICase1[i, j]
deriv(ICase1[, 2:8]) <- (gamma_E_t * E2[i, j] * prob_hosp_t_mult[i, j]) - (gamma_ICase_t * ICase1[i, j]) - gamma_vaccine_t[j] * ICase1[i, j] + gamma_vaccine_t[j - 1] * ICase1[i, j - 1]

deriv(ICase2[, 1]) <- (gamma_ICase_t * ICase1[i, j]) - (gamma_ICase_t * ICase2[i, j]) - gamma_vaccine_t[j] * ICase2[i, j]
deriv(ICase2[, 2:8]) <- (gamma_ICase_t * ICase1[i, j]) - (gamma_ICase_t * ICase2[i, j]) - gamma_vaccine_t[j] * ICase2[i, j] + gamma_vaccine_t[j - 1] * ICase2[i, j - 1]

output(ICase[]) <- sum(ICase1[i, ]) + sum(ICase2[i, ])
dim(ICase) <- 17
################################################################################

### IOxGetLive (IOxGetLive1 & IOxGetLive2): Get oxygen, go on to survive #######
dim(IOxGetLive1) <- c(17, 8)
dim(IOxGetLive2) <- c(17, 8)

IOxGetLive1_0[, ] <- user()
dim(IOxGetLive1_0) <- c(17, 8)
initial(IOxGetLive1[, ]) <- IOxGetLive1_0[i, j]

IOxGetLive2_0[, ] <- user()
dim(IOxGetLive2_0) <- c(17, 8)
initial(IOxGetLive2[, ]) <- IOxGetLive2_0[i, j]

gamma_get_ox_survive[] <- user() # rate of progression through requiring oxygen compartment conditional on getting oxygen and surviving
dim(gamma_get_ox_survive) <- user()
tt_dur_get_ox_survive[] <- user()
dim(tt_dur_get_ox_survive) <- length(gamma_get_ox_survive)
gamma_get_ox_survive_t <- interpolate(tt_dur_get_ox_survive, gamma_get_ox_survive, "constant")


deriv(IOxGetLive1[, 1]) <- (gamma_ICase_t * ICase2[i, j] * (1 - prob_severe_multi[i]) * p_oxygen * (1 - prob_non_severe_death_treatment[i])) - (gamma_get_ox_survive_t * IOxGetLive1[i, j]) - gamma_vaccine_t[j] * IOxGetLive1[i, j]
deriv(IOxGetLive1[, 2:8]) <- (gamma_ICase_t * ICase2[i, j] * (1 - prob_severe_multi[i]) * p_oxygen * (1 - prob_non_severe_death_treatment[i])) - (gamma_get_ox_survive_t * IOxGetLive1[i, j]) - gamma_vaccine_t[j] * IOxGetLive1[i, j] + gamma_vaccine_t[j - 1] * IOxGetLive1[i, j - 1]

deriv(IOxGetLive2[, 1]) <- (gamma_get_ox_survive_t * IOxGetLive1[i, j]) - (gamma_get_ox_survive_t * IOxGetLive2[i, j]) - gamma_vaccine_t[j] * IOxGetLive2[i, j]
deriv(IOxGetLive2[, 2:8]) <- (gamma_get_ox_survive_t * IOxGetLive1[i, j]) - (gamma_get_ox_survive_t * IOxGetLive2[i, j]) - gamma_vaccine_t[j] * IOxGetLive2[i, j] + gamma_vaccine_t[j - 1] * IOxGetLive2[i, j - 1]

################################################################################

### IOxGetDie (IOxGetDie1 & IOxGetDie2): Get oxygen go on to die ###############
dim(IOxGetDie1) <- c(17, 8)
dim(IOxGetDie2) <- c(17, 8)

IOxGetDie1_0[, ] <- user()
dim(IOxGetDie1_0) <- c(17, 8)
initial(IOxGetDie1[, ]) <- IOxGetDie1_0[i, j]

IOxGetDie2_0[, ] <- user()
dim(IOxGetDie2_0) <- c(17, 8)
initial(IOxGetDie2[, ]) <- IOxGetDie2_0[i, j]

gamma_get_ox_die[] <- user() # rate of progression through requiring oxygen compartment conditional on getting oxygen and dying
dim(gamma_get_ox_die) <- user()
tt_dur_get_ox_die[] <- user()
dim(tt_dur_get_ox_die) <- length(gamma_get_ox_die)
gamma_get_ox_die_t <- interpolate(tt_dur_get_ox_die, gamma_get_ox_die, "constant")

deriv(IOxGetDie1[, ]) <- (gamma_ICase_t * ICase2[i, j] * (1 - prob_severe_multi[i]) * p_oxygen * prob_non_severe_death_treatment[i]) - gamma_get_ox_die_t * IOxGetDie1[i, j]

deriv(IOxGetDie2[, ]) <-  (gamma_get_ox_die_t * IOxGetDie1[i, j]) - (gamma_get_ox_die_t * IOxGetDie2[i, j])

################################################################################

### IOxNotGetLive (IOxNotGetLive1 & IOxNotGetLive2): Do not get oxygen, go on to survive #######
dim(IOxNotGetLive1) <- c(17, 8)
dim(IOxNotGetLive2) <- c(17, 8)

IOxNotGetLive1_0[, ] <- user()
dim(IOxNotGetLive1_0) <- c(17, 8)
initial(IOxNotGetLive1[, ]) <- IOxNotGetLive1_0[i, j]

IOxNotGetLive2_0[, ] <- user()
dim(IOxNotGetLive2_0) <- c(17, 8)
initial(IOxNotGetLive2[, ]) <- IOxNotGetLive2_0[i, j]

gamma_not_get_ox_survive <- user() # rate of progression through requiring oxygen compartment conditional on not getting oxygen and surviving

deriv(IOxNotGetLive1[, 1]) <-  (gamma_ICase_t * ICase2[i, j] * (1 - prob_severe_multi[i]) * (1 - p_oxygen) * (1 - prob_non_severe_death_no_treatment[i])) - (gamma_not_get_ox_survive * IOxNotGetLive1[i, j]) - gamma_vaccine_t[j] * IOxNotGetLive1[i, j]
deriv(IOxNotGetLive1[, 2:8]) <- (gamma_ICase_t * ICase2[i, j] * (1 - prob_severe_multi[i]) * (1 - p_oxygen) * (1 - prob_non_severe_death_no_treatment[i])) - (gamma_not_get_ox_survive * IOxNotGetLive1[i, j]) - gamma_vaccine_t[j] * IOxNotGetLive1[i, j] + gamma_vaccine_t[j - 1] * IOxNotGetLive1[i, j - 1]

deriv(IOxNotGetLive2[, 1]) <- (gamma_not_get_ox_survive * IOxNotGetLive1[i, j]) - (gamma_not_get_ox_survive * IOxNotGetLive2[i, j]) - gamma_vaccine_t[j] * IOxNotGetLive2[i, j]
deriv(IOxNotGetLive2[, 2:8]) <- (gamma_not_get_ox_survive * IOxNotGetLive1[i, j]) - (gamma_not_get_ox_survive * IOxNotGetLive2[i, j]) - gamma_vaccine_t[j] * IOxNotGetLive2[i, j] + gamma_vaccine_t[j - 1] * IOxNotGetLive2[i, j - 1]

################################################################################

### IOxNotGetDie (IOxNotGetDie1 & IOxNotGetDie2): Do not get oxygen, go on to die #######
dim(IOxNotGetDie1) <- c(17, 8)
dim(IOxNotGetDie2) <- c(17, 8)

IOxNotGetDie1_0[, ] <- user()
dim(IOxNotGetDie1_0) <- c(17, 8)
initial(IOxNotGetDie1[, ]) <- IOxNotGetDie1_0[i, j]

IOxNotGetDie2_0[, ] <- user()
dim(IOxNotGetDie2_0) <- c(17, 8)
initial(IOxNotGetDie2[, ]) <- IOxNotGetDie2_0[i, j]

gamma_not_get_ox_die <- user() # rate of progression through requiring oxygen compartment conditional on not getting oxygen and dying

deriv(IOxNotGetDie1[, ]) <- (gamma_ICase_t * ICase2[i, j] * (1 - prob_severe_multi[i]) * (1 - p_oxygen) * prob_non_severe_death_no_treatment[i]) - (gamma_not_get_ox_die * IOxNotGetDie1[i, j])

deriv(IOxNotGetDie2[, ]) <- (gamma_not_get_ox_die * IOxNotGetDie1[i, j]) - (gamma_not_get_ox_die * IOxNotGetDie2[i, j])

################################################################################

### IMVGetLive (IMVGetLive1 & IMVGetLive2): Get mechanical ventilation, go on to live ########
dim(IMVGetLive1) <- c(17, 8)
dim(IMVGetLive2) <- c(17, 8)

IMVGetLive1_0[, ] <- user()
dim(IMVGetLive1_0) <- c(17, 8)
initial(IMVGetLive1[, ]) <- IMVGetLive1_0[i, j]

IMVGetLive2_0[, ] <- user()
dim(IMVGetLive2_0) <- c(17, 8)
initial(IMVGetLive2[, ]) <- IMVGetLive2_0[i, j]

gamma_get_mv_survive[] <- user() # rate of progression through requiring mechanical ventilation compartment conditional on getting ventilation and surviving
dim(gamma_get_mv_survive) <- user()
tt_dur_get_mv_survive[] <- user()
dim(tt_dur_get_mv_survive) <- length(gamma_get_mv_survive)
gamma_get_mv_survive_t <- interpolate(tt_dur_get_mv_survive, gamma_get_mv_survive, "constant")

deriv(IMVGetLive1[, 1]) <- (gamma_ICase_t * ICase2[i, j] * prob_severe_multi[i] * p_ventilation * (1 - prob_severe_death_treatment[i])) - (gamma_get_mv_survive_t * IMVGetLive1[i, j]) - gamma_vaccine_t[j] * IMVGetLive1[i, j]
deriv(IMVGetLive1[, 2:8]) <- (gamma_ICase_t * ICase2[i, j] * prob_severe_multi[i] * p_ventilation * (1 - prob_severe_death_treatment[i])) - (gamma_get_mv_survive_t * IMVGetLive1[i, j]) - gamma_vaccine_t[j] * IMVGetLive1[i, j] + gamma_vaccine_t[j - 1] * IMVGetLive1[i, j - 1]

deriv(IMVGetLive2[, 1]) <- (gamma_get_mv_survive_t * IMVGetLive1[i, j]) - (gamma_get_mv_survive_t * IMVGetLive2[i, j]) - gamma_vaccine_t[j] * IMVGetLive2[i, j]
deriv(IMVGetLive2[, 2:8]) <- (gamma_get_mv_survive_t * IMVGetLive1[i, j]) - (gamma_get_mv_survive_t * IMVGetLive2[i, j]) - gamma_vaccine_t[j] * IMVGetLive2[i, j] + gamma_vaccine_t[j - 1] * IMVGetLive2[i, j - 1]

################################################################################

### IMVGetDie (IMVGetDie1 & IMVGetDie2): Get mechanical ventilation, go on to die ########
dim(IMVGetDie1) <- c(17, 8)
dim(IMVGetDie2) <- c(17, 8)

IMVGetDie1_0[, ] <- user()
dim(IMVGetDie1_0) <- c(17, 8)
initial(IMVGetDie1[, ]) <- IMVGetDie1_0[i, j]

IMVGetDie2_0[, ] <- user()
dim(IMVGetDie2_0) <- c(17, 8)
initial(IMVGetDie2[, ]) <- IMVGetDie2_0[i, j]

gamma_get_mv_die[] <- user() # rate of progression through requiring mechanical ventilation compartment conditional on getting ventilation and dying
dim(gamma_get_mv_die) <- user()
tt_dur_get_mv_die[] <- user()
dim(tt_dur_get_mv_die) <- length(gamma_get_mv_die)
gamma_get_mv_die_t <- interpolate(tt_dur_get_mv_die, gamma_get_mv_die, "constant")

deriv(IMVGetDie1[, ]) <- (gamma_ICase_t * ICase2[i, j] * prob_severe_multi[i] * p_ventilation * prob_severe_death_treatment[i]) - (gamma_get_mv_die_t * IMVGetDie1[i, j])

deriv(IMVGetDie2[, ]) <- (gamma_get_mv_die_t * IMVGetDie1[i, j]) - (gamma_get_mv_die_t * IMVGetDie2[i, j])

################################################################################

### IMVNotGetLive (IMVNotGetLive1 & IMVNotGetLive2): Do no get mechanical ventilation, go on to live ########
dim(IMVNotGetLive1) <- c(17, 8)
dim(IMVNotGetLive2) <- c(17, 8)

IMVNotGetLive1_0[, ] <- user()
dim(IMVNotGetLive1_0) <- c(17, 8)
initial(IMVNotGetLive1[, ]) <- IMVNotGetLive1_0[i, j]

IMVNotGetLive2_0[, ] <- user()
dim(IMVNotGetLive2_0) <- c(17, 8)
initial(IMVNotGetLive2[, ]) <- IMVNotGetLive2_0[i, j]

gamma_not_get_mv_survive <- user() # rate of progression through requiring mechanical ventilation compartment conditional on not getting ventilation and surviving

deriv(IMVNotGetLive1[, 1]) <- (gamma_ICase_t * ICase2[i, j] * prob_severe_multi[i] * (1 - p_ventilation) * (1 - prob_severe_death_no_treatment[i])) - (gamma_not_get_mv_survive * IMVNotGetLive1[i, j]) - gamma_vaccine_t[j] * IMVNotGetLive1[i, j]
deriv(IMVNotGetLive1[, 2:8]) <- (gamma_ICase_t * ICase2[i, j] * prob_severe_multi[i] * (1 - p_ventilation) * (1 - prob_severe_death_no_treatment[i])) - (gamma_not_get_mv_survive * IMVNotGetLive1[i, j]) - gamma_vaccine_t[j] * IMVNotGetLive1[i, j] + gamma_vaccine_t[j - 1] * IMVNotGetLive1[i, j - 1]

deriv(IMVNotGetLive2[, 1]) <- (gamma_not_get_mv_survive * IMVNotGetLive1[i, j]) - (gamma_not_get_mv_survive * IMVNotGetLive2[i, j]) - gamma_vaccine_t[j] * IMVNotGetLive2[i, j]
deriv(IMVNotGetLive2[, 2:8]) <- (gamma_not_get_mv_survive * IMVNotGetLive1[i, j]) - (gamma_not_get_mv_survive * IMVNotGetLive2[i, j]) - gamma_vaccine_t[j] * IMVNotGetLive2[i, j] + gamma_vaccine_t[j - 1] * IMVNotGetLive2[i, j - 1]

################################################################################

### IMVNotGetDie (IMVNotGetDie1 & IMVNotGetDie2): Do no get mechanical ventilation, go on to die ########
dim(IMVNotGetDie1) <- c(17, 8)
dim(IMVNotGetDie2) <- c(17, 8)

IMVNotGetDie1_0[, ] <- user()
dim(IMVNotGetDie1_0) <- c(17, 8)
initial(IMVNotGetDie1[, ]) <- IMVNotGetDie1_0[i, j]

IMVNotGetDie2_0[, ] <- user()
dim(IMVNotGetDie2_0) <- c(17, 8)
initial(IMVNotGetDie2[, ]) <- IMVNotGetDie2_0[i, j]

gamma_not_get_mv_die <- user() # rate of progression through requiring mechanical ventilation compartment conditional on not getting ventilation and dying


deriv(IMVNotGetDie1[, ]) <- (gamma_ICase_t * ICase2[i, j] * prob_severe_multi[i] * (1 - p_ventilation) * prob_severe_death_no_treatment[i]) - (gamma_not_get_mv_die * IMVNotGetDie1[i, j])

deriv(IMVNotGetDie2[, ]) <- (gamma_not_get_mv_die * IMVNotGetDie1[i, j]) - (gamma_not_get_mv_die * IMVNotGetDie2[i, j])

################################################################################

### IRec (IRec1 & IRec2): Recovering from ICU ##################################
dim(IRec1) <- c(17, 8)
dim(IRec2) <- c(17, 8)
dim(IRec) <- 17

IRec1_0[, ] <- user()
dim(IRec1_0) <- c(17, 8)
initial(IRec1[, ]) <- IRec1_0[i, j]

IRec2_0[, ] <- user()
dim(IRec2_0) <- c(17, 8)
initial(IRec2[, ]) <- IRec2_0[i, j]

gamma_rec <- user() # rate of progression through post-ICU recovery compartment

deriv(IRec1[, 1]) <- (gamma_get_mv_survive_t * IMVGetLive2[i, j]) - (gamma_rec * IRec1[i, j]) - gamma_vaccine_t[j] * IRec1[i, j]
deriv(IRec1[, 2:8]) <- (gamma_get_mv_survive_t * IMVGetLive2[i, j]) - (gamma_rec * IRec1[i, j]) - gamma_vaccine_t[j] * IRec1[i, j] + gamma_vaccine_t[j - 1] * IRec1[i, j - 1]

deriv(IRec2[, 1]) <- (gamma_rec * IRec1[i, j]) - (gamma_rec * IRec2[i, j]) - gamma_vaccine_t[j] * IRec2[i, j]
deriv(IRec2[, 2:8]) <- (gamma_rec * IRec1[i, j]) - (gamma_rec * IRec2[i, j]) - gamma_vaccine_t[j] * IRec2[i, j] + gamma_vaccine_t[j - 1] * IRec2[i, j - 1]

output(IRec[]) <- sum(IRec1[i, ]) + sum(IRec2[i, ])
################################################################################

### D: Dead ####################################################################
dim(D) <- c(17, 8)

D_0[, ] <- user()
dim(D_0) <- c(17, 8)
initial(D[, ]) <- D_0[i, j]

deriv(D[, 1:8]) <- (gamma_get_ox_die_t * IOxGetDie2[i, j]) + (gamma_not_get_ox_die * IOxNotGetDie2[i, j]) + (gamma_get_mv_die_t * IMVGetDie2[i, j]) + (gamma_not_get_mv_die * IMVNotGetDie2[i, j])
################################################################################

################################################################################
### Vaccination capacity #######################################################
################################################################################
# Vaccination
# Vaccine prioritisation coverage matrix
N_prioritisation_steps <- user()
vaccine_coverage_mat[, ] <- user()
dim(vaccine_coverage_mat) <- c(N_prioritisation_steps, 17)

# Generating Vaccine Efficacy Over Time
vaccine_efficacy_infection_t[, ] <- interpolate(tt_vaccine_efficacy_infection, vaccine_efficacy_infection, "constant")
dim(vaccine_efficacy_infection_t) <- c(17, 8)
tt_vaccine_efficacy_infection[] <- user()
vaccine_efficacy_infection[, , ] <- user()
dim(tt_vaccine_efficacy_infection) <- user()
dim(vaccine_efficacy_infection) <- c(length(tt_vaccine_efficacy_infection), 17, 8)

gamma_vaccine[,] <- user() # Vector of vaccine progression parameters by vaccination status (only effects rate of waning)
tt_dur_vaccine[] <- user()
dim(gamma_vaccine) <- c(length(tt_dur_vaccine), 8)
dim(tt_dur_vaccine) <- user()
gamma_vaccine_t[] <- interpolate(tt_dur_vaccine, gamma_vaccine, "constant")
dim(gamma_vaccine_t) <- 8

# Interpolation of vaccination rate over time
### Original
# t_primary_doses <- interpolate(tt_primary_doses, primary_doses, "constant")
# tt_primary_doses[] <- user()
# primary_doses[] <- user()
# dim(tt_primary_doses) <- user()
# dim(primary_doses) <- length(tt_primary_doses)

### Original
# nuisance_dose_length_parameter <- length(primary_doses) + 1
# second_dose_delay <- user()
# t_second_doses <- interpolate(tt_second_doses, second_doses, "constant")
# tt_second_doses[2:nuisance_dose_length_parameter] <- tt_primary_doses[i - 1] + second_dose_delay
# tt_second_doses[1] <- 0
# second_doses[2:nuisance_dose_length_parameter] <- primary_doses[i - 1]
# second_doses[1] <- 0
# dim(second_doses) <- length(tt_second_doses)
# dim(tt_second_doses) <- nuisance_dose_length_parameter

### Original
# t_booster_doses <- interpolate(tt_booster_doses, booster_doses, "constant")
# tt_booster_doses[] <- user()
# booster_doses[] <- user()
# dim(tt_booster_doses) <- user()
# dim(booster_doses) <- length(tt_booster_doses)


### Altered so that we can have BPSV and Spec have different second dose delays (even though they're both put in the model as the same vaccine series for diff age groups)
runtime <- user()

primary_doses[] <- user()     # Calculate outside (including the time-variable delay to protection)
dim(primary_doses) <- runtime

second_doses[] <- user()      # Calculate outside (including the second_dose_delay)
dim(second_doses) <- runtime

booster_doses[] <- user()       # Calculate outside (including the time-variable delay to protection)
dim(booster_doses) <- user()

######

vaccine_booster_follow_up_coverage[] <- user()
dim(vaccine_booster_follow_up_coverage) <- 17

vaccine_booster_initial_coverage[] <- user()
dim(vaccine_booster_initial_coverage) <- 17

# Track the number who have received each type of vaccine in each age group, maybe make this just those who can be vaccinated?
dose_pops[, ] <- sum(S[i, j:8]) + sum(E1[i, j:8]) + sum(E2[i, j:8]) + sum(IMild[i, j:8]) + sum(ICase1[i, j:8]) + sum(ICase2[i, j:8]) +
  sum(IMVGetLive1[i, j:8]) + sum(IMVGetLive2[i, j:8]) +
  sum(IMVGetDie1[i, j:8]) + sum(IMVGetDie2[i, j:8]) + sum(IMVNotGetLive1[i, j:8]) + sum(IMVNotGetLive2[i, j:8]) + sum(IMVNotGetDie1[i, j:8]) + sum(IMVNotGetDie2[i, j:8]) +
  sum(IOxGetLive1[i, j:8]) + sum(IOxGetLive2[i, j:8]) + sum(IOxGetDie1[i, j:8]) + sum(IOxGetDie2[i, j:8]) + sum(IOxNotGetLive1[i, j:8]) + sum(IOxNotGetLive2[i, j:8]) +
  sum(IOxNotGetDie1[i, j:8]) + sum(IOxNotGetDie2[i, j:8]) +
  sum(IRec1[i, j:8]) + sum(IRec2[i, j:8]) +
  sum(R1[i, j:8]) + sum(R2[i, j:8])
dim(dose_pops) <- c(17, 7) # 1 is everyone alive, 2 is all with first dose, 3 is all with second dose, 4,5 is all initial waned, 6 is all with first booster dose, 7 is all who've waned

# number of people at each vaccination level
# useful to see vaccination levels, also used in outputs
vaccination_cov[,] <- dose_pops[i, j] - dose_pops[i, j + 1]
dim(vaccination_cov) <- c(17, 6)
#1: Unvaccinated, 2: First Dose, 3: Second Dose, 4,5:Waned Second Dose, 6:Boosted

# Calculate priorisation step for the first doses (THIS BREAKS IF NOT INCLUSIVE OF PREVIOUS TARGETS, WRITE A CHECK!!!)
target_met_matrix[, ] <- (vaccine_coverage_mat[i, j] * dose_pops[j, 1]) <= (dose_pops[j, 2] + 1)
dim(target_met_matrix) <- c(N_prioritisation_steps, 17)
target_met_column[] <- sum(target_met_matrix[i, ]) == 17
dim(target_met_column) <- c(N_prioritisation_steps)
prioritisation_step <- if (sum(target_met_column) < N_prioritisation_steps) sum(target_met_column) + 1 else N_prioritisation_steps

# Calculate number of people available to vaccinate for first/second/boosters
target_pop_first[] <- max(((vaccine_coverage_mat[as.integer(prioritisation_step), i] * dose_pops[i, 1]) - dose_pops[i, 2]), 0)
dim(target_pop_first) <- 17
# number of doses
primary_first[] <- min(primary_doses[as.integer(t)] * target_pop_first[i] / max(sum(target_pop_first) * (vaccination_cov[i,1]), 1), 1)
dim(primary_first) <- 17

## Original
primary_second[] <- min(second_doses[as.integer(t)] / max(sum(vaccination_cov[i, 2]), 1), 1)
dim(primary_second) <- 17

## Consider changing to make age specific??
# primary_second[] <- min(second_doses[as.integer(t)] / max(sum(vaccination_cov[, 2]), 1), 1)
# dim(primary_second) <- 17

#first boosters:
#ISSUE: I'd like to stack them up but Odins dependency check sucks, maybe mention to Rich
vaccination_cov_eligible[, ] <- vaccination_cov[i, j + 2] * vaccine_booster_initial_coverage[i]
dim(vaccination_cov_eligible) <- c(17, 3)
eligible_for_first_booster[] <- sum(vaccination_cov_eligible[i, ])
dim(eligible_for_first_booster) <- 17
booster_first[] <- min(booster_doses[as.integer(t)] / max(sum(eligible_for_first_booster[]), 1), 1) * vaccine_booster_initial_coverage[i]
dim(booster_first) <- 17

a_initial_boosted[] <- booster_first[i] * eligible_for_first_booster[i]
dim(a_initial_boosted) <- 17
remaining_boosters <- booster_doses[as.integer(t)] - sum(a_initial_boosted)
eligible_for_follow_up_booster[] <- dose_pops[i, 7] * vaccine_booster_follow_up_coverage[i]
dim(eligible_for_follow_up_booster) <- 17
booster_second[] <- min(remaining_boosters / max(sum(eligible_for_follow_up_booster[]), 1), 1) * vaccine_booster_follow_up_coverage[i]
#we can just multiply because these are flat percentages really
dim(booster_second) <- 17


################################################################################
################################################################################

################################################################################
### Hospital and ICU capacity ##################################################
################################################################################
## Interpolation for Hospital and ICU Capacity
hosp_bed_capacity <- interpolate(tt_hosp_beds, hosp_beds, "constant")
tt_hosp_beds[] <- user()
hosp_beds[] <- user()
dim(tt_hosp_beds) <- user()
dim(hosp_beds) <- length(tt_hosp_beds)

ICU_bed_capacity <- interpolate(tt_ICU_beds, ICU_beds, "constant")
tt_ICU_beds[] <- user()
ICU_beds[] <- user()
dim(tt_ICU_beds) <- user()
dim(ICU_beds) <- length(tt_ICU_beds)

# Generating prob_hosp Over Time
prob_hosp_t[, ] <- interpolate(tt_vaccine_efficacy_disease, prob_hosp, "constant")
prob_hosp_t_mult[, ] <- prob_hosp_multiplier_t * prob_hosp_t[i, j]
dim(prob_hosp_t) <- c(17, 8) # probability of requiring hospitalisation by age and vaccination status at time t
dim(prob_hosp_t_mult) <- c(17, 8) # probability of requiring hospitalisation by age and vaccination status at time t
tt_vaccine_efficacy_disease[] <- user()
prob_hosp[, , ] <- user()
dim(tt_vaccine_efficacy_disease) <- user()
dim(prob_hosp) <- c(length(tt_vaccine_efficacy_disease), 17, 8)

# Interpolation for prob_hosp_multiplier
prob_hosp_multiplier_t <- interpolate(tt_prob_hosp_multiplier, prob_hosp_multiplier, "constant")
tt_prob_hosp_multiplier[] <- user()
dim(tt_prob_hosp_multiplier) <- user()
dim(prob_hosp_multiplier) <- length(tt_prob_hosp_multiplier)
prob_hosp_multiplier[] <- user() # rate of progression through recovered compartment (loss of naturally acquired immunity)

# Probability of severe symptoms with multiplier
prob_severe[] <- user() # probability of severe disease (requiring mechanical ventilation) by age
dim(prob_severe) <- 17

prob_severe_multiplier[] <- user() # flat modifer to severity
dim(prob_severe_multiplier) <- length(tt_prob_severe_multiplier)
tt_prob_severe_multiplier[] <- user()
dim(tt_prob_severe_multiplier) <- user()

# interpolate severity
prob_severe_multiplier_t <- interpolate(tt_prob_severe_multiplier, prob_severe_multiplier, "constant")
# calculate the new severity
prob_severe_multi[] <- prob_severe_multiplier_t * prob_severe[i]
dim(prob_severe_multi) <- 17

prob_non_severe_death_treatment[] <- user() # probability of dying from non-severe disease (i.e. requiring oxygen but not mechanical ventilation) by age given you receive appropriate treatment (proxy here is whether a general hospital bed is available)
dim(prob_non_severe_death_treatment) <- 17

prob_non_severe_death_no_treatment[] <- user() # probability of dying from non-severe disease (i.e. requiring oxygen but not mechanical ventilation) by age given you do NOT receive appropriate treatment (proxy here is whether a general hospital bed is available)
dim(prob_non_severe_death_no_treatment) <- 17

prob_severe_death_treatment[] <- user() # probability of dying from severe disease (i.e. requiring mechanical ventilation) by age given you receive appropriate treatment (proxy here is whether an ICU bed is available)
dim(prob_severe_death_treatment) <- 17

prob_severe_death_no_treatment[] <- user() # probability of dying from severe disease (i.e. requiring mechanical ventilation) by age given you do NOT receive appropriate treatment (proxy here is whether an ICU bed is available)
dim(prob_severe_death_no_treatment) <- 17

rel_infectiousness[] <- user() # Relative infectiousness of age categories relative to maximum infectiousness age category
dim(rel_infectiousness) <- 17

rel_infectiousness_vaccinated[, ] <- user() # Relative infectiousness of vaccinated age categories relative to maximum infectiousness age category
dim(rel_infectiousness_vaccinated) <- c(17, 8)

# Infections Requiring Oxygen (a general Hosptial Bed)
hosp_occ <- sum(IOxGetLive1) + sum(IOxGetLive2) - gamma_get_ox_survive_t * sum(IOxGetLive2) + sum(IOxGetDie1) + sum(IOxGetDie2) - gamma_get_ox_die_t * sum(IOxGetDie2) + sum(IRec1) + sum(IRec2) - gamma_rec * sum(IRec2) # Summing number of infections in compartments that use general hospital beds
number_requiring_Ox[, ] <- gamma_ICase_t * ICase2[i, j] * (1 - prob_severe_multi[i])
dim(number_requiring_Ox) <- c(17, 8)
total_number_requiring_ox <- sum(number_requiring_Ox)

p_oxygen <- if ((total_number_requiring_ox <= (hosp_bed_capacity - hosp_occ)) || total_number_requiring_ox <= 0) 1 else (hosp_bed_capacity - hosp_occ) / total_number_requiring_ox

# Infections Requiring Mechanical Ventilation (an ICU Bed)
ICU_occ <- sum(IMVGetLive1) + sum(IMVGetLive2) - gamma_get_mv_survive_t * sum(IMVGetLive2) + sum(IMVGetDie1) + sum(IMVGetDie2) - gamma_get_mv_die_t * sum(IMVGetDie2) # Summing number of infections in compartments that use ICU beds
number_requiring_IMV[, ] <- gamma_ICase_t * ICase2[i, j] * prob_severe_multi[i]
dim(number_requiring_IMV) <- c(17, 8)
total_number_requiring_IMV <- sum(number_requiring_IMV)

p_ventilation <- if (total_number_requiring_IMV <= (ICU_bed_capacity - ICU_occ) || total_number_requiring_IMV <= 0) 1 else (ICU_bed_capacity - ICU_occ) / total_number_requiring_IMV

################################################################################
################################################################################

################################################################################
### FOI and contact matrix #####################################################
################################################################################
# Generating Force of Infection
m[, ] <- interpolate(tt_matrix, mix_mat_set, "constant")
dim(m) <- c(17, 17)
tt_matrix[] <- user()
mix_mat_set[, , ] <- user()
dim(tt_matrix) <- user()
dim(mix_mat_set) <- c(length(tt_matrix), 17, 17)

# Interpolation for beta
beta <- interpolate(tt_beta, beta_set, "constant")
tt_beta[] <- user()
beta_set[] <- user()
dim(tt_beta) <- user()
dim(beta_set) <- length(tt_beta)

# Generating Force of Infection
temp_rel[, ] <- (IMild[i, j] * rel_infectiousness_vaccinated[i, j]) + (ICase1[i, j] * rel_infectiousness_vaccinated[i, j]) + (ICase2[i, j] * rel_infectiousness_vaccinated[i, j])
temp[] <- sum(temp_rel[i, ])
dim(temp_rel) <- c(17, 8)
dim(temp) <- c(17)

s_ij[, ] <- m[i, j] * temp[j] * rel_infectiousness[j]
dim(s_ij) <- c(17, 17)

lambda[] <- beta * sum(s_ij[i, ])
dim(lambda) <- 17
################################################################################
################################################################################

################################################################################
### Output #####################################################################
################################################################################
# Hospital occupancy and demand
output(hospital_occupancy[]) <- sum(IOxGetLive1[i, ]) + sum(IOxGetLive2[i, ]) + sum(IOxGetDie1[i, ]) + sum(IOxGetDie2[i, ]) + sum(IRec1[i, ]) + sum(IRec2[i, ])
dim(hospital_occupancy) <- 17

output(ICU_occupancy[]) <- sum(IMVGetLive1[i, ]) + sum(IMVGetLive2[i, ]) + sum(IMVGetDie1[i, ]) + sum(IMVGetDie2[i, ])
dim(ICU_occupancy) <- 17

output(hospital_demand[]) <- sum(IOxGetLive1[i, ]) + sum(IOxGetLive2[i, ]) + sum(IOxGetDie1[i, ]) + sum(IOxGetDie2[i, ]) + sum(IRec1[i, ]) + sum(IRec2[i, ]) + sum(IOxNotGetLive1[i, ]) + sum(IOxNotGetLive2[i, ]) + sum(IOxNotGetDie1[i, ]) + sum(IOxNotGetDie2[i, ])
dim(hospital_demand) <- 17

output(ICU_demand[]) <- sum(IMVGetLive1[i, ]) + sum(IMVGetLive2[i, ]) + sum(IMVGetDie1[i, ]) + sum(IMVGetDie2[i, ]) + sum(IMVNotGetLive1[i, ]) + sum(IMVNotGetLive2[i, ]) + sum(IMVNotGetDie1[i, ]) + sum(IMVNotGetDie2[i, ])
dim(ICU_demand) <- 17

# Number in hospital or ICU compartments
output(IICU[]) <- sum(IMVGetLive1[i, ]) + sum(IMVGetLive2[i, ]) + sum(IMVGetDie1[i, ]) + sum(IMVGetDie2[i, ]) + sum(IMVNotGetLive1[i, ]) + sum(IMVNotGetLive2[i, ]) + sum(IMVNotGetDie1[i, ]) + sum(IMVNotGetDie2[i, ])
dim(IICU) <- 17

output(IHospital[]) <- sum(IOxGetLive1[i, ]) + sum(IOxGetLive2[i, ]) + sum(IOxGetDie1[i, ]) + sum(IOxGetDie2[i, ]) + sum(IOxNotGetLive1[i, ]) + sum(IOxNotGetLive2[i, ]) + sum(IOxNotGetDie1[i, ]) + sum(IOxNotGetDie2[i, ])
dim(IHospital) <- 17

# Hospitalisations
deriv(hospitalisations_cumu[, ]) <- number_requiring_IMV[i, j] * p_ventilation + number_requiring_Ox[i, j] * p_oxygen
dim(hospitalisations_cumu) <- c(17, 8)
initial(hospitalisations_cumu[, ]) <- 0

# Hospitalisation demand
deriv(hospitalisation_demand_cumu[, ]) <- number_requiring_IMV[i, j] + number_requiring_Ox[i, j]
dim(hospitalisation_demand_cumu) <- c(17, 8)
initial(hospitalisation_demand_cumu[, ]) <- 0

# Deaths
output(deaths_cumu[, ]) <- D[i, j]
dim(deaths_cumu) <- c(17, 8)

# Infections
deriv(infections_cumu[, ]) <- (lambda[i] * vaccine_efficacy_infection_t[i, j] * S[i, j])
dim(infections_cumu) <- c(17, 8)
initial(infections_cumu[, ]) <- 0

output(first_doses_given[]) <- primary_first[i] * (S[i, 1] + E1[i, 1] + E2[i, 1] + R1[i, 1] + R2[i, 1])
dim(first_doses_given) <- 17

output(second_doses_given[]) <- primary_second[i] * (S[i, 2] + E1[i, 2] + E2[i, 2] + R1[i, 2] + R2[i, 2])
dim(second_doses_given) <- 17

output(booster_doses_given[]) <- booster_first[i] * (sum(S[i, 3:5]) + sum(E1[i, 3:5]) +
                                                       sum(E2[i, 3:5]) + sum(R1[i, 3:5]) + sum(R2[i, 3:5])) +
  booster_second[i] * (sum(S[i, 7:8]) + sum(E1[i, 7:8]) +
                         sum(E2[i, 7:8]) + sum(R1[i, 7:8]) + sum(R2[i, 7:8]))
dim(booster_doses_given) <- 17


# Unvaccinated
output(unvaccinated[]) <- vaccination_cov[i, 1]
dim(unvaccinated) <- 17
# Vaccinated First Dose
output(vaccinated_first_dose[]) <- dose_pops[i, 2]
dim(vaccinated_first_dose) <- 17
# Vaccinated Second Dose
output(vaccinated_second_dose[]) <- dose_pops[i, 3]
dim(vaccinated_second_dose) <- 17
# Vaccinated Second Dose Waned
output(vaccinated_second_waned[]) <-  dose_pops[i, 4] -  dose_pops[i, 6]
dim(vaccinated_second_waned) <- 17
# Vaccinated booster Dose
output(vaccinated_booster_dose[]) <- dose_pops[i, 6]
dim(vaccinated_booster_dose) <- 17
# Waned
output(vaccinated_booster_waned[]) <- dose_pops[i, 7]
dim(vaccinated_booster_waned) <- 17

output(N[]) <- dose_pops[i, 1] + sum(D[i, ])
dim(N) <- 17
################################################################################
