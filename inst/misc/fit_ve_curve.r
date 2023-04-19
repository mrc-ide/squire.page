library(tidyverse)
library(odin)
library(ggpubr)
library(dfoptim)
library(forcats)
set.seed(1000101)
#assumed numbers, booster restores efficacy NOTE Ideally we'd find better numbers, might be too many categories
master_ves <- tribble(
  ~dose, ~endpoint, ~efficacy,
  "First", "Infection", 0.55,
  "Second", "Infection", 0.75,
  "Booster", "Infection", 0.80,
  "First", "Hospitalisation", 0.75,
  "Second", "Hospitalisation", 0.90,
  "Booster", "Hospitalisation", 0.95
)

ab_params <- list(
  ni50 = log10(0.052),
  ns50 = log10(0.01),
  k = 2.8,
  hl_s = 33,
  hl_l = 580,
  period_s = 86
)

##Fit Waning Curves
simulate_time <- 3*365
calc_eff_gen <- odin({
  initial(C1) <- 1
  initial(C2) <- 0
  initial(C3) <- 0
  deriv(C1) <- -w_1*C1
  deriv(C2) <- w_1*C1 - w_2*C2
  deriv(C3) <- w_2*C2
  output(ve_d) <- (C1 * ved + C2 * ved_2 + C3 * ved_3)
  output(ve_i) <- (C1 * vei + C2 * vei_2 + C3 * vei_3)
  w_1 <- user()
  w_2 <- user()
  ved <- user()
  ved_2 <- user()
  ved_3 <- user()
  vei <- user()
  vei_2 <- user()
  vei_3 <- user()
})

simulate_ab <- function(t, initial_ab, h_s, h_l, t_s) {
  pi1 <- -log(2)/h_s
  pi2 <- -log(2)/h_l
  initial_ab * (
    (exp(pi1 * t + pi2 * t_s) + exp(pi1 * t_s + pi2 * t))/
    (exp(pi1 * t_s) + exp(pi2 * t_s))
  )
}

ab_to_ve <- function(ab, n50, k){
  1/(1 + exp(-k * (log10(ab) - n50)))
}

err_lines <- function(l1, l2){
  #sum((l1 - l2)^2/l1)
  sum(((l1 - l2)/l1)^2)
  #scale it so the lower values have more weight
}

calculate_ve <- function(p1, p2, p3) {
  ved <- p1
  ved_2 <- ved * p2
  ved_3 <- ved_2 * p3
  c(ved, ved_2, ved_3)
}

first_doses <- master_ves %>%
  filter(dose == "First") %>%
  group_by(dose) %>%
  group_split() %>%
  map(function(df){
    parameter_infection <- df %>% filter(endpoint == "Infection") %>% pull(efficacy)
    parameter_hospitalisation <- df %>% filter(endpoint == "Hospitalisation") %>% pull(efficacy)
    dose <- unique(df$dose)
    tibble(
      parameter = c("pV_i", "pV_d"),
      value = c(parameter_infection, (parameter_hospitalisation - parameter_infection)/(1 - parameter_infection)),
      dose = dose
    )
  })

other_doses <- master_ves %>%
  filter(dose != "First") %>%
  group_by(dose) %>%
  group_split() %>%
  map(function(df){
    parameter_infection <- df %>% filter(endpoint == "Infection") %>% pull(efficacy)
    parameter_hospitalisation <- df %>% filter(endpoint == "Hospitalisation") %>% pull(efficacy)

    dose <- unique(df$dose)

    message(paste0(dose))

    t_plot <- seq(0, simulate_time, length.out = 100)

    calc_eff <- calc_eff_gen$new(user = list(
      ved = parameter_hospitalisation,
      vei = parameter_infection,
      ved_2 = 1,
      ved_3 = 1,
      vei_2 = 1,
      vei_3 = 1,
      w_1 = 0,
      w_2 = 0
    ))

    #need to calculate initial dose level
    #assume ve is 30 days after dose
    t_measure <- 30
    #just assume the initia AB is
    err_func <- function(initial_ab) {
      ab_t_measure <- simulate_ab(t_measure, initial_ab, ab_params$hl_s, ab_params$hl_l, ab_params$period_s)
      ve_i <- ab_to_ve(ab_t_measure, ab_params$ni50, ab_params$k)
      ve_d <- ab_to_ve(ab_t_measure, ab_params$ns50, ab_params$k)
      err_lines(parameter_infection, ve_i) +
        err_lines(parameter_hospitalisation, ve_d)
    }
    res <- optimize(err_func, interval = c(0, 10), maximum = FALSE)
    initial_ab <- res$minimum
    abs <- simulate_ab(0:simulate_time, initial_ab,  ab_params$hl_s,  ab_params$hl_l,  ab_params$period_s)
    ve_i <- ab_to_ve(abs, ab_params$ni50, ab_params$k)
    ve_d <- ab_to_ve(abs, ab_params$ns50, ab_params$k)
    #plot for diagnostics
    initial_ab_plot <- tibble(
      t = 0:(2*t_measure),
      Hospitalisation = ve_d[1:(2*t_measure + 1)],
      Infection = ve_i[1:(2*t_measure + 1)]
    ) %>% 
      pivot_longer(cols = -t, names_to = "Endpoint", values_to = "Efficacy") %>%
      ggplot(aes(x = t, y = Efficacy, colour = Endpoint)) + 
      geom_line(show.legend = FALSE) +
      geom_point(data = tibble(
        t = t_measure,
        Efficacy = c(parameter_infection, parameter_hospitalisation),
        Endpoint = c("Infection", "Hospitalisation")
      ), shape = "x", size = 5, show.legend = FALSE) +
      ggpubr::theme_pubclean() +
      labs(
        title = "Fit of initial Immune Response level, X shows the input efficacies",
        y = "Vaccine Efficacy", x = "Days Since Dose"
      )

    #scale for break through infection
    ve_d <- (ve_d - ve_i)/(1 - ve_i)

    err_func <- function(pars) {
      veds <- calculate_ve(pars[3], pars[4], pars[5])
      veis <- calculate_ve(pars[6], pars[7], pars[8])
      calc_eff$set_user(
        user = list(
          w_1 = pars[1],
          w_2 = pars[2],
          ved = veds[1],
          ved_2 = veds[2],
          ved_3 = veds[3],
          vei = veis[1],
          vei_2 = veis[2],
          vei_3 = veis[3]
        )
      )
      mod_value <- calc_eff$run(t = seq(0, simulate_time))
      log(sqrt(err_lines(ve_d, mod_value[, "ve_d"]) + err_lines(ve_i, mod_value[, "ve_i"])))
    }
    lower = list(
      w_1 = 1/(3*365),
      w_2 = 1/(3*365),
      ved_1 = 0,
      ved_2 = 0,
      ved_3 = 0,
      vei_1 = 0,
      vei_2 = 0,
      vei_3 = 0
    )
    upper = list(
      w_1 = 1/30,
      w_2 = 1/30,
      ved_1 = 1,
      ved_2 = 1,
      ved_3 = 1,
      vei_1 = 1,
      vei_2 = 1,
      vei_3 = 1
    )
    par =  list(
      w_1 = 1/365,
      w_2 = 1/365,
      ved_1 = 0.5,
      ved_2 = 0.5,
      ved_3 = 0.5,
      vei_1 = 0.5,
      vei_2 = 0.5,
      vei_3 = 0.5
    )
    res <- dfoptim::nmkb(unlist(par), fn = err_func, lower = unlist(lower), upper = unlist(upper), control = list(maxfeval = 5000))
    if(res$convergence != 0){
      stop(res$message)
    }

    #add randomness (should do it in fitting really)
    out <- map(1, function(i){
      pars <- res$par
      veds <- calculate_ve(pars[3], pars[4], pars[5])
      veis <- calculate_ve(pars[6], pars[7], pars[8])
      tibble(
        value = c(pars[1:2], veds, veis),
        parameter = names(par)
      )
    })[[1]]

    #plot
    t_plot <- 0:simulate_time
    calc_eff$set_user(
        user = list(
            w_1 = out$value[1],
            w_2 = out$value[2],
            ved = out$value[3],
            ved_2 = out$value[4],
            ved_3 = out$value[5],
            vei = out$value[6],
            vei_2 = out$value[7],
            vei_3 = out$value[8]
        )
    )
    mod_value <- calc_eff$run(t = t_plot)
    p <- tibble(
        t = rep(t_plot, 2),
        value = c(mod_value[, "ve_d"], mod_value[, "ve_i"]),
        endpoint = c(rep("Hospitalisation", length(t_plot)), rep("Infection", length(t_plot)))
    ) %>%
    mutate(
        model = "Booster Model"
    ) %>%
      rbind(
        tibble(
          t = t_plot,
          value = ve_d,
          endpoint = "Hospitalisation",
          model = "AB Process"
        )
      ) %>%
      rbind(
        tibble(
          t = t_plot,
          value = ve_i,
          endpoint = "Infection",
          model = "AB Process"
        )
      ) %>%
      ggplot(aes(x = t, y = value, color = endpoint, linetype = model)) +
        geom_line() +
        labs(y = "Vaccine Efficacy", x = "Days Since Dose", title = paste0("Dose: ", dose), linetype = "Model", colour = "Endpoint") +
        ggpubr::theme_pubclean() + scale_alpha(guide = 'none') +
        ylim(c(0, 1))

    out$dose <- dose

    list(
      out = out,
      plot = p,
      initial_ab_plot = initial_ab_plot
    )
})
#split into plots and data
plots <- map(other_doses, ~.x$plot)
initial_ab_plot <- map(other_doses, ~.x$initial_ab_plot)
other_doses <- map(other_doses, ~.x$out)

##Calibration plots
p <- ggarrange(
    ggarrange(plotlist = plots, ncol = 2, common.legend = TRUE),
    ggarrange(plotlist = initial_ab_plot, ncol = 2),
    ncol = 1, heights = c(0.6, 0.2)
)

pdf("vignettes/fitted_value_plot.pdf", width = 15, height = 15)
print(p)
dev.off()

rm(p, plots, initial_ab_plot)

first_doses <- map_dfr(first_doses, ~.x)
other_doses <- map_dfr(other_doses, ~.x)

other_doses <- other_doses %>%
  mutate(
    parameter = case_when(
      parameter == "w_1" & dose == "Second" ~ "fw_1",
      parameter == "w_2" & dose == "Second"~ "fw_2",
      parameter == "ved_1" & dose == "Second" ~ "fV_d_1",
      parameter == "ved_2" & dose == "Second" ~ "fV_d_2",
      parameter == "ved_3" & dose == "Second" ~ "fV_d_3",
      parameter == "vei_1" & dose == "Second" ~ "fV_i_1",
      parameter == "vei_2" & dose == "Second" ~ "fV_i_2",
      parameter == "vei_3" & dose == "Second" ~ "fV_i_3",
      parameter == "w_1" & dose == "Booster" ~ "bw_1",
      parameter == "w_2" & dose == "Booster"~ "bw_2",
      parameter == "ved_1" & dose == "Booster" ~ "bV_d_1",
      parameter == "ved_2" & dose == "Booster" ~ "bV_d_2",
      parameter == "ved_3" & dose == "Booster" ~ "bV_d_3",
      parameter == "vei_1" & dose == "Booster" ~ "bV_i_1",
      parameter == "vei_2" & dose == "Booster" ~ "bV_i_2",
      parameter == "vei_3" & dose == "Booster" ~ "bV_i_3",
      TRUE ~ parameter
    )
  )

#combine
out <- first_doses %>%
  rbind(other_doses)

efficacies <- out %>%
    transmute(
        Endpoint = case_when(
            str_detect(parameter, "_i") ~ "Protection against Infection",
            str_detect(parameter, "_d") ~ "Protection against Hospitalisation",
            TRUE ~ NA
        ),
        Compartment = str_remove(parameter, "_i|_d"),
        value = value,
        Compartment = fct_relevel(Compartment, "pV", "fV_1", "fV_2", "fV_3", "bV_1", "bV_2", "bV_3")
    ) %>%
    filter(!is.na(Endpoint)) %>% 
    pivot_wider(names_from = Endpoint, values_from = value) %>% 
    arrange(Compartment)

waning_rates <- out %>% 
    filter(str_detect(parameter, "w_")) %>% 
    pull(value, parameter)
print(efficacies)
print(waning_rates)
print(1/waning_rates)

