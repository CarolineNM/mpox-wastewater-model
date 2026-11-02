# ==================================================
# Script: 04_run_counterfactuals.R
# Purpose: Run counterfactual scenario analyses
# Project: Mpox_clean
# Updated: 2026-04-22
# ==================================================

rm(list = ls())
set.seed(123)

# ---- Packages ----
library(tidyverse)
library(coda)
library(readxl)
library(patchwork)


# ---- Transit time posterior summary ----
# min = 0.4, median = 3.22, mode (~Q3) = 4.69, upper = 7.26


# ---- Population structure and contact rates ----
dat <- data.frame(
  setting      = "Vancouver",
  Sex_groups   = c("0-4", "4-22", "22-100"),
  mean_rate    = c(1.31, 9.13, 39.84),
  mean_rate_dy = c(1.31, 9.13, 39.84) / 180,
  Prop         = c(0.60, 0.35, 0.05),
  Pop          = c(0.60, 0.35, 0.05) * 26100
)

# convert 6-month sexual contact rates to daily rates
contact <- dat$mean_rate / 180
contact_dy <- dat$mean_rate_dy

# normalized daily contact rates for distributing exposure risk
norm_contact_dy <- contact_dy / sum(contact_dy)

# population size by activity group
N <- c(15660, 9135, 1305)

# simulation horizon
burn_in_timesteps <- 30
T <- 169 + burn_in_timesteps


# ---- Wastewater flow data ----
ww_flow <- read_excel("data/raw/case_data_V2.xlsx", sheet = "flowrate")

flow_all <- ww_flow %>% select(sumflow)
flow_alldaily <- as.numeric(unlist(flow_all))

# prepend burn-in period and convert litres to mL
flow_mlalldaily <- c(
  rep(0, burn_in_timesteps),
  flow_alldaily * 1e3
)


# ---- Case data ----
case_dat <- read_excel("data/raw/case_data_V2.xlsx", sheet = "cases")

cases_obsb <- case_dat$total_cases
Cumulat_inc <- case_dat$Cumulat_Inc


# ---- Reporting delay data ----
rdelay <- read_excel("data/raw/case_data_V2.xlsx", sheet = "rdelay")
summary(rdelay)


# ---- Wastewater observations ----
ww_dat <- read_excel("data/raw/case_data_V2.xlsx", sheet = "dailyWW")

ww_std <- ww_dat %>% select(log10_cp_per_person_per_day)
ww_obs <- as.numeric(unlist(ww_std$log10_cp_per_person_per_day))

# align observation days with simulation timeline
obs_days <- ww_dat$Day + burn_in_timesteps


# ---- Transit time scenarios ----
# posterior-anchored scenarios

transit_times <- c(
  fast     = 1.96,   # lower scenario
  moderate = 4.89,   # upper-middle scenario
  slow     = 7.35    # upper scenario
)

# posterior median
baseline_transit <- 3.54

# multiplicative factors relative to baseline
transit_factors <- transit_times / baseline_transit


# ---- Scenario wrapper ----
run_transit_scenarios <- function(
  params, init, contact, N, T,
  factors = transit_factors
) {
  
  out <- list()
  
  for (f in factors) {
    
    p2 <- params
    p2$transit_mean <- params$transit_mean * f
    
    name <- if (abs(f - 1) < 1e-9) {
      "Baseline transit"
    } else {
      paste0("Longer transit: +", round((f - 1) * 100), "%")
    }
    
    out[[name]] <- mysimulate_model(
      p2, init, contact, N, T,
      vax_scenario = "baseline"
    )
  }
  
  out
}

# ---- Asymptomatic infection scenarios ----
# Multiplicative increases in the asymptomatic proportion (m)
m_factors <- c(1.20, 1.60, 1.80)

run_asym_scenarios <- function(
  params, init_vals, contact, N, T,
  factors
) {
  
  results <- vector("list", length(factors))
  
  for (i in seq_along(factors)) {
    
    p_mod <- params
    
    # scale asymptomatic proportion
    p_mod$m <- p_mod$m * factors[i]
    
    out <- mysimulate_model(
      p_mod, init_vals, contact, N, T,
      vax_scenario = "baseline"
    )
    
    out$factor <- factors[i]
    results[[i]] <- out
  }
  
  # scenario labels, e.g. "Asymptomatic: +20%"
  labels <- paste0(
    "Asymptomatic: ",
    ifelse(factors >= 1, "+", ""),
    round((factors - 1) * 100),
    "%"
  )
  
  setNames(results, labels)
}


# ---- Reporting delay scenarios ----
# Alternative delays (days) between infection detection and reporting
reporting_delays <- c(5, 10, 14)

run_reporting_delay_scenarios <- function(
  params, init, contact, N, T,
  delays = c(5, 10, 14)
) {
  
  out <- list()
  
  for (d in delays) {
    
    p_mod <- params
    
    # override baseline reporting delay
    p_mod$report_lag <- d
    
    name <- paste0("Reporting delay: ", d, " days")
    
    out[[name]] <- mysimulate_model(
      p_mod, init, contact, N, T,
      vax_scenario = "baseline"
    )
  }
  
  out
}


# ---- Convert posterior draw to simulation parameters ----
# Input: one posterior sample row
extract_params <- function(draw_row) {
  
  list(
    
    # transmission and mixing
    beta = draw_row["beta"],
    kappa = draw_row["kappa"],
    phi = draw_row["phi"],
    
    # vaccine effectiveness
    Vea = draw_row["Vea"],
    rsa = 1 - draw_row["Vea"],
    Veb = draw_row["Veb"],
    rsb = 1 - draw_row["Veb"],
    
    # reporting process
    report_frac = draw_row["report_frac"],
    report_lag = 7,
    
    # asymptomatic proportion
    m = draw_row["m"],
    mva = draw_row["m"],
    mvb = draw_row["m"],
    
    # progression rates
    delta = 1 / draw_row["delta_inv"],
    epsilon = 1 / 2,
    theta = 1 / (draw_row["theta_invall"] / 3),
    omega = 1 / (draw_row["omega_invall"] / 3),
    
    deltava = 1 / draw_row["delta_inv"],
    epsilonva = 1 / 2,
    thetava = 1 / (draw_row["theta_invall"] / 3),
    omegava = 1 / (draw_row["omega_invall"] / 3),
    
    deltavb = 1 / draw_row["delta_inv"],
    epsilonvb = 1 / 2,
    thetavb = 1 / (draw_row["theta_invall"] / 3),
    omegavb = 1 / (draw_row["omega_invall"] / 3),
    
    # wastewater and shedding parameters
    mu = 0.18,
    tau = 0.776,
    
    shed_I1 = 10^12.3,
    shed_I2 = 10^11.5,
    shed_I3 = 10^10.2,
    
    shed_A1 = 10^12.3 * 0.776,
    shed_A2 = 10^11.5 * 0.776,
    shed_A3 = 10^10.2 * 0.776,
    
    shed_alpha = 10^12.3 * 0.776,
    
    shed_P1a = 10^12.3 * 0.776,
    shed_P1b = 10^12.3 * 0.776,
    
    shed_A1a = 10^12.3 * 0.776,
    shed_A2a = 10^11.5 * 0.776,
    shed_A3a = 10^10.2 * 0.776,
    
    shed_A1b = 10^12.3 * 0.776,
    shed_A2b = 10^11.5 * 0.776,
    shed_A3b = 10^10.2 * 0.776,
    
    shed_I1a = 10^12.3,
    shed_I2a = 10^11.5,
    shed_I3a = 10^10.2,
    
    shed_I1b = 10^12.3,
    shed_I2b = 10^11.5,
    shed_I3b = 10^10.2,
    
    # transport and observation
    transit_mean = draw_row["transit_time_mean"],
    transit_cv = draw_row["transit_time_cv"],
    tau_ww = draw_row["tau_ww"],
    mult = exp(draw_row["log_mult"]),
    
    # vaccination delays
    va_delay = 14,
    vb_delay = 14,
    
    # wastewater system inputs
    flow_mlalldaily = flow_mlalldaily,
    wwtp_population = 1278020
  )
}

# ---- Load fitted model output ----
# Use posterior samples from the combined model fit
load(file.path("outputs", "fits", "combined_fit.RData"))
mcmc_matrixallcom <- as.matrix(do.call(rbind, Combined_finaltestfbnovurbb))
# select one posterior draw for simulation
one_draw <- mcmc_matrixallcom[1000, ]
params <- extract_params(one_draw)

# ---- Simulation model ----
# Forward simulation function with flexible vaccination scenarios
mysimulate_model <- function(params, init, contact, N, T, 
                             vax_scenario = "baseline", 
                             v1_frac = 1, 
                             v2_frac = 1) {
  G <- length(N)
  
  alpha_v <- rep(0, T)
  mu_v <- rep(0, T)
  
  if (vax_scenario == "baseline") {
    alpha_v <- init$alpha_v
    mu_v <- init$mu_v
    
    # Force vaccinated compartments to shed like unvaccinated
    params$shed_P1a <- params$shed_alpha
    params$shed_P1b <- params$shed_alpha
    
    params$shed_I1a <- params$shed_I1
    params$shed_I2a <- params$shed_I2
    params$shed_I3a <- params$shed_I3
    params$shed_I1b <- params$shed_I1
    params$shed_I2b <- params$shed_I2
    params$shed_I3b <- params$shed_I3
    
    params$shed_A1a <- params$shed_A1
    params$shed_A2a <- params$shed_A2
    params$shed_A3a <- params$shed_A3
    params$shed_A1b <- params$shed_A1
    params$shed_A2b <- params$shed_A2
    params$shed_A3b <- params$shed_A3
    
  } else if (vax_scenario == "reduced_shedding") {
    alpha_v <- init$alpha_v
    mu_v    <- init$mu_v
    
    # Apply dose-specific shedding multipliers (remaining fraction)
    # Dose 1
    params$shed_P1a <- params$shed_alpha * v1_frac
    params$shed_I1a <- params$shed_I1    * v1_frac
    params$shed_I2a <- params$shed_I2    * v1_frac
    params$shed_I3a <- params$shed_I3    * v1_frac
    params$shed_A1a <- params$shed_A1    * v1_frac
    params$shed_A2a <- params$shed_A2    * v1_frac
    params$shed_A3a <- params$shed_A3    * v1_frac
    
    # Dose 2
    params$shed_P1b <- params$shed_alpha * v2_frac
    params$shed_I1b <- params$shed_I1    * v2_frac
    params$shed_I2b <- params$shed_I2    * v2_frac
    params$shed_I3b <- params$shed_I3    * v2_frac
    params$shed_A1b <- params$shed_A1    * v2_frac
    params$shed_A2b <- params$shed_A2    * v2_frac
    params$shed_A3b <- params$shed_A3    * v2_frac
    
    message(sprintf("[Scenario: reduced shedding] v1_frac=%.2f, v2_frac=%.2f", 
                    v1_frac, v2_frac))
    
  } else {
    stop("vax_scenario must be 'baseline' or 'reduced_shedding'")
  }
  
  
  # Extract parameters
  va_delay <- params$va_delay
  vb_delay <- params$vb_delay
  beta <- params$beta
  kappa <- params$kappa
  phi <- params$phi
  rsa <- params$rsa
  rsb <- params$rsb
  report_frac <- params$report_frac
  report_lag <- params$report_lag
  m <- params$m
  mva <- params$mva
  mvb <- params$mvb
  delta <- params$delta
  epsilon <- params$epsilon
  theta <- params$theta
  omega <- params$omega
  deltava <- params$deltava
  epsilonva <- params$epsilonva
  thetava <- params$thetava
  omegava <- params$omegava
  deltavb <- params$deltavb
  epsilonvb <- params$epsilonvb
  thetavb <- params$thetavb
  omegavb <- params$omegavb
  mu <- params$mu
  tau <- params$tau
  shed_I1 <- params$shed_I1
  shed_I2 <- params$shed_I2
  shed_I3 <- params$shed_I3
  shed_A1 <- params$shed_A1
  shed_A2 <- params$shed_A2
  shed_A3 <- params$shed_A3
  shed_alpha <- params$shed_alpha
  shed_P1a<-params$shed_P1a
  shed_P1b<-params$shed_P1b
  shed_A1a <- params$shed_A1a
  shed_A2a <- params$shed_A2a
  shed_A3a <- params$shed_A3a
  shed_A1b <- params$shed_A1b
  shed_A2b <- params$shed_A2b
  shed_A3b <- params$shed_A3b
  shed_I1a <- params$shed_I1a
  shed_I2a <- params$shed_I2a
  shed_I3a <- params$shed_I3a
  shed_I1b <- params$shed_I1b
  shed_I2b <- params$shed_I2b
  shed_I3b <- params$shed_I3b
  transit_mean <- params$transit_mean
  transit_cv <- params$transit_cv
  tau_ww <- params$tau_ww
  mult <- params$mult
  flow <- params$flow_mlalldaily
  wwtp_population <- params$wwtp_population
  
  
  va_queue <- array(0, dim = c(G, va_delay, T))
  vb_queue <- array(0, dim = c(G, vb_delay, T))
  
  # For vaccine effectiveness kick-in
  Va_eff <- matrix(0, G, T)
  Vb_eff <- matrix(0, G, T)
  shed_P<-matrix(0, G, T)
  shed_A<-matrix(0, G, T)
  shed_I<-matrix(0, G, T)
  P_total<-matrix(0, G, T)
  A_total<-matrix(0, G, T)
  I_total<-matrix(0, G, T)
  # Compartment initialization function
  comp_matrix <- function() matrix(0, G, T)
  S <- comp_matrix(); E <- comp_matrix(); P <- comp_matrix()
  A1 <- comp_matrix(); A2 <- comp_matrix(); A3 <- comp_matrix()
  I1 <- comp_matrix(); I2 <- comp_matrix(); I3 <- comp_matrix()
  R <- comp_matrix(); Cuminc <- comp_matrix()
  Va <- comp_matrix(); Eva <- comp_matrix(); Pva <- comp_matrix()
  A1va <- comp_matrix(); A2va <- comp_matrix(); A3va <- comp_matrix()
  I1va <- comp_matrix(); I2va <- comp_matrix(); I3va <- comp_matrix()
  Vb <- comp_matrix(); Evb <- comp_matrix(); Pvb <- comp_matrix()
  A1vb <- comp_matrix(); A2vb <- comp_matrix(); A3vb <- comp_matrix()
  I1vb <- comp_matrix(); I2vb <- comp_matrix(); I3vb <- comp_matrix()
  prevalence <- comp_matrix(); lambda_det <- comp_matrix()
  
  new_reported_cases <- matrix(0, G, T)
  total_lambda <- numeric(T)
  active_infected <- numeric(T)
  total_new_cases <- numeric(T)
  total_Cuminc <- numeric(T)
  daily_shedding<-numeric(T)
  
  # Mixing matrix setup
  mix_ass <- diag(1, G)
  mix_prp <- outer(contact, contact * N) / sum(contact * N)
  mix <- array(0, dim = c(G, G, T))
  
  for (g in 1:G) {
    for (g_ in 1:G) {
      mix[g, g_, 1] <- kappa * mix_ass[g, g_] + (1 - kappa) * mix_prp[g, g_]
    }
  }
  
  # Initial states
  S[, 1] <- init$S
  E[, 1] <- init$E
  P[, 1] <- init$P
  A1[, 1] <- init$A1
  I1[, 1] <- init$I1
  
  # Calculate initial prevalence and lambda
  for (g in 1:G) {
    prevalence[g, 1] <- sum(c(P[g,1], A1[g,1], I1[g,1])) / N[g]
    lambda_det[g, 1] <- sum(beta * contact[g] * mix[g, , 1] * prevalence[, 1])
  }
  total_lambda[1] <- sum(lambda_det[, 1])
  
  # ---- Initialize transition matrices for each group ----
  new_Es     <- matrix(0, G, T)
  new_Ps     <- matrix(0, G, T)
  new_A1s    <- matrix(0, G, T)
  new_A2s    <- matrix(0, G, T)
  new_A3s    <- matrix(0, G, T)
  new_RAs    <- matrix(0, G, T)
  new_I1s    <- matrix(0, G, T)
  new_I2s    <- matrix(0, G, T)
  new_I3s    <- matrix(0, G, T)
  new_RIs    <- matrix(0, G, T)
  
  new_Vas    <- matrix(0, G, T)
  new_Evas   <- matrix(0, G, T)
  new_Pvas   <- matrix(0, G, T)
  new_A1vas  <- matrix(0, G, T)
  new_A2vas  <- matrix(0, G, T)
  new_A3vas  <- matrix(0, G, T)
  new_RAvas  <- matrix(0, G, T)
  new_I1vas  <- matrix(0, G, T)
  new_I2vas  <- matrix(0, G, T)
  new_I3vas  <- matrix(0, G, T)
  new_RIvas  <- matrix(0, G, T)
  
  new_Vbs    <- matrix(0, G, T)
  new_Evbs   <- matrix(0, G, T)
  new_Pvbs   <- matrix(0, G, T)
  new_A1vbs  <- matrix(0, G, T)
  new_A2vbs  <- matrix(0, G, T)
  new_A3vbs  <- matrix(0, G, T)
  new_RAvbs  <- matrix(0, G, T)
  new_I1vbs  <- matrix(0, G, T)
  new_I2vbs  <- matrix(0, G, T)
  new_I3vbs  <- matrix(0, G, T)
  new_RIvbs  <- matrix(0, G, T)
  
  
  # Dynamics of the model (deterministic-like)
  for (t in 2:T) {
    
    # Final mixing matrix for three groups only
    for (g in 1:3) {
      for (g_ in 1:3) {
        mix[g, g_, t] <- kappa * mix_ass[g, g_] + (1 - kappa) * mix_prp[g, g_]
      }
    }
    
    # Calculate prevalence once per time step (deterministic)
    # Calculate prevalence once per time step (deterministic)
    for (g in 1:3) {
      prevalence[g, t] <- ((tau*P[g, t-1]) + I1[g, t-1] + I2[g, t-1] + I3[g, t-1]+
                             (tau*A1[g, t-1]) + (tau*A2[g, t-1]) + (tau*A3[g, t-1])+
                             (tau*Pva[g, t-1])+I1va[g, t-1] + I2va[g, t-1] + I3va[g, t-1]+
                             (tau*A1va[g, t-1])+ (tau*A2va[g, t-1])+ (tau*A3va[g, t-1])+
                             (tau*Pvb[g, t-1])+I1vb[g, t-1] + I2vb[g, t-1] + I3vb[g, t-1]+
                             (tau*A1vb[g, t-1])+ (tau*A2vb[g, t-1])+ (tau*A3vb[g, t-1]))/ N[g]
      
    }
    
    # Calculate force of infection (lambda_det) for each group
    for (g in 1:3) {
      lambda_det[g, t] <- sum(beta * contact[g] * mix[g, 1:3, t] * prevalence[1:3, t])
    }
    total_lambda[t] <- sum(lambda_det[1:3, t])
    
    # Deterministic transitions instead of binomial draws
    for (g in 1:3) {
      new_Es[g, t-1] <- lambda_det[g, t] * S[g, t-1]     # Deterministic S to E transition
      new_Ps[g, t-1] <- delta * E[g, t-1]             # Deterministic E to P transition
      new_A1s[g, t-1] <- epsilon * m * P[g, t-1]      # Deterministic P to A1
      new_A2s[g, t-1] <- omega * A1[g, t-1]           # Deterministic A1 to A2
      new_A3s[g, t-1] <- omega * A2[g, t-1]           # Deterministic A2 to A3
      new_RAs[g, t-1] <- omega * A3[g, t-1]           # Deterministic A3 to R
      new_I1s[g, t-1] <- epsilon * (1 - m) * P[g, t-1] # Deterministic P to I1
      new_I2s[g, t-1] <- theta * I1[g, t-1]           # Deterministic I1 to I2
      new_I3s[g, t-1] <- theta * I2[g, t-1]           # Deterministic I2 to I3
      new_RIs[g, t-1] <- theta * I3[g, t-1]           # Deterministic I3 to R
      
      
      
      # Deterministic transitions for 1st dose vaccination
      
      new_Vas[g, t-1] <- alpha_v[t] * S[g, t-1]           # Deterministic S to Va transition.alpha is vaccination rate
      # Shift the Sv_queue (buffer of people waiting for protection)
      va_queue[g, 1, t] <- new_Vas[g, t-1]
      for (d in 2:va_delay) {
        va_queue[g, d, t] <- va_queue[g, d-1, t-1]
      }
      
      
      Va_eff[g, t] <- va_queue[g, va_delay, t-1]       # Individuals whose protection kicks in today
      new_Evas[g, t-1] <- lambda_det[g, t] *rsa* Va[g, t-1]     # Deterministic Va to Eva transition.Lambda*1-vaccine effectiveness(ve)
      new_Pvas[g, t-1] <- deltava * Eva[g, t-1]             # Deterministic Eva to Pva transition
      new_A1vas[g, t-1] <- epsilonva * mva * Pva[g, t-1]      # Deterministic Pva to A1va
      new_A2vas[g, t-1] <- omegava * A1va[g, t-1]           # Deterministic A1va to A2va
      new_A3vas[g, t-1] <- omegava * A2va[g, t-1]           # Deterministic A2va to A3va
      new_RAvas[g, t-1] <- omegava * A3va[g, t-1]           # Deterministic A3va to Rva
      new_I1vas[g, t-1] <- epsilonva * (1 - mva) * Pva[g, t-1] # Deterministic Pva to I1va
      new_I2vas[g, t-1] <- thetava * I1va[g, t-1]           # Deterministic I1va to I2va
      new_I3vas[g, t-1] <- thetava * I2va[g, t-1]           # Deterministic I2va to I3va
      new_RIvas[g, t-1] <- thetava * I3va[g, t-1]           # Deterministic I3va to Rv
      
      # Deterministic transitions for 2nd dose vaccination
      
      new_Vbs[g, t-1] <- mu_v[t] * Va[g, t-1]              # Deterministic Va to Vb transition.mu is vaccination rate
      # Shift the Sv_queue (buffer of people waiting for protection)
      vb_queue[g, 1, t] <- new_Vbs[g, t-1]
      for (d in 2:vb_delay) {
        vb_queue[g, d, t] <- vb_queue[g, d-1, t-1]
      }
      
      Vb_eff[g, t] <- vb_queue[g, vb_delay, t-1]
      new_Evbs[g, t-1] <- lambda_det[g, t] *rsb* Vb[g, t-1]     # Deterministic Vb to Evb transition.Lambda*1-vaccine effectiveness(ve)
      new_Pvbs[g, t-1] <- deltavb * Evb[g, t-1]             # Deterministic Evb to Pvb transition
      new_A1vbs[g, t-1] <- epsilonvb * mvb * Pvb[g, t-1]      # Deterministic Pvb to A1vb
      new_A2vbs[g, t-1] <- omegavb * A1vb[g, t-1]           # Deterministic A1vb to A2vb
      new_A3vbs[g, t-1] <- omegavb * A2vb[g, t-1]           # Deterministic A2vb to A3vb
      new_RAvbs[g, t-1] <- omegavb * A3vb[g, t-1]           # Deterministic A3vb to RAvb
      new_I1vbs[g, t-1] <- epsilonvb * (1 - mvb) * Pvb[g, t-1] # Deterministic Pvb to I1vb
      new_I2vbs[g, t-1] <- thetavb * I1vb[g, t-1]           # Deterministic I1vb to I2vb
      new_I3vbs[g, t-1] <- thetavb * I2vb[g, t-1]           # Deterministic I2vb to I3vb
      new_RIvbs[g, t-1] <- thetavb * I3vb[g, t-1]           # Deterministic I3vb to RIvb
      
      
      
      # Update compartments deterministically at each time step
      
      S[g, t] <- max(0, S[g, t-1] - new_Es[g, t-1]-new_Vas[g, t-1])
      E[g, t] <- max(0, E[g, t-1] + new_Es[g, t-1] - new_Ps[g, t-1])
      P[g, t] <- max(0, P[g, t-1] + new_Ps[g, t-1] - new_A1s[g, t-1] - new_I1s[g, t-1])
      A1[g, t] <- max(0, A1[g, t-1] + new_A1s[g, t-1] - new_A2s[g, t-1])
      A2[g, t] <- max(0, A2[g, t-1] + new_A2s[g, t-1] - new_A3s[g, t-1])
      A3[g, t] <- max(0, A3[g, t-1] + new_A3s[g, t-1] - new_RAs[g, t-1])
      I1[g, t] <- max(0, I1[g, t-1] + new_I1s[g, t-1] - new_I2s[g, t-1])
      I2[g, t] <- max(0, I2[g, t-1] + new_I2s[g, t-1] - new_I3s[g, t-1])
      I3[g, t] <- max(0, I3[g, t-1] + new_I3s[g, t-1] - new_RIs[g, t-1])
      R[g, t] <- max(0, R[g, t-1] + new_RAs[g, t-1] + new_RIs[g, t-1]+new_RAvas[g, t-1]+new_RIvas[g, t-1]+
                       new_RAvbs[g, t-1]+new_RIvbs[g, t-1])
      
      Cuminc[g, t] <- Cuminc[g, t-1] + new_Es[g, t-1]+ new_Evas[g, t-1] + new_Evbs[g, t-1]
      
      # Add comps with 1st dose vaccination Va
      
      Va[g, t] <- max(0, Va[g, t-1] + Va_eff[g, t] - new_Evas[g, t-1] - new_Vbs[g, t-1])
      Eva[g, t] <- max(0, Eva[g, t-1] + new_Evas[g, t-1] - new_Pvas[g, t-1])
      Pva[g, t] <- max(0, Pva[g, t-1] + new_Pvas[g, t-1] - new_A1vas[g, t-1] - new_I1vas[g, t-1])
      A1va[g, t] <- max(0, A1va[g, t-1] + new_A1vas[g, t-1] - new_A2vas[g, t-1])
      A2va[g, t] <- max(0, A2va[g, t-1] + new_A2vas[g, t-1] - new_A3vas[g, t-1])
      A3va[g, t] <- max(0, A3va[g, t-1] + new_A3vas[g, t-1] - new_RAvas[g, t-1])
      I1va[g, t] <- max(0, I1va[g, t-1] + new_I1vas[g, t-1] - new_I2vas[g, t-1])
      I2va[g, t] <- max(0, I2va[g, t-1] + new_I2vas[g, t-1] - new_I3vas[g, t-1])
      I3va[g, t] <- max(0, I3va[g, t-1] + new_I3vas[g, t-1] - new_RIvas[g, t-1])
      
      
      # Add comps with 2nd dose vaccination Vb
      
      Vb[g, t] <- max(0, Vb[g, t-1] + Vb_eff[g, t]-new_Evbs[g, t-1])
      Evb[g, t] <- max(0, Evb[g, t-1] + new_Evbs[g, t-1] - new_Pvbs[g, t-1])
      Pvb[g, t] <- max(0, Pvb[g, t-1] + new_Pvbs[g, t-1] - new_A1vbs[g, t-1] - new_I1vbs[g, t-1])
      A1vb[g, t] <- max(0, A1vb[g, t-1] + new_A1vbs[g, t-1] - new_A2vbs[g, t-1])
      A2vb[g, t] <- max(0, A2vb[g, t-1] + new_A2vbs[g, t-1] - new_A3vbs[g, t-1])
      A3vb[g, t] <- max(0, A3vb[g, t-1] + new_A3vbs[g, t-1] - new_RAvbs[g, t-1])
      I1vb[g, t] <- max(0, I1vb[g, t-1] + new_I1vbs[g, t-1] - new_I2vbs[g, t-1])
      I2vb[g, t] <- max(0, I2vb[g, t-1] + new_I2vbs[g, t-1] - new_I3vbs[g, t-1])
      I3vb[g, t] <- max(0, I3vb[g, t-1] + new_I3vbs[g, t-1] - new_RIvbs[g, t-1])
      
      
      # Calculate reported cases based on all infectious compartments with reporting lag
      # Calculate reported cases based on new_I1s and new_I1vas only
      step <- function(x) {
        ifelse(x >= 1, 1, 0)
      }
      new_reported_cases[g, t] <- step(t - report_lag) * report_frac *
        (new_I1s[g, max(1, t - report_lag)]+new_I1vas[g, max(1, t - report_lag)]+
           new_I1vbs[g, max(1, t - report_lag)])
      
    }
    
    # Total new reported cases at time t(incidence at t)
    total_new_cases[t] <- sum(new_reported_cases[1:3, t])
    
    total_Cuminc[t] <- sum(Cuminc[1:3, t])
    
    # Update active infections
    active_infected[t] <-
      sum(E[1:3, t]) +
      sum(P[1:3, t]) +
      sum(A1[1:3, t]) + sum(A2[1:3, t]) + sum(A3[1:3, t]) +
      sum(I1[1:3, t]) + sum(I2[1:3, t]) + sum(I3[1:3, t]) +
      sum(Pva[1:3, t]) + sum(Pvb[1:3, t]) +
      sum(A1va[1:3, t]) + sum(A2va[1:3, t]) + sum(A3va[1:3, t]) +
      sum(I1va[1:3, t]) + sum(I2va[1:3, t]) + sum(I3va[1:3, t]) +
      sum(A1vb[1:3, t]) + sum(A2vb[1:3, t]) + sum(A3vb[1:3, t]) +
      sum(I1vb[1:3, t]) + sum(I2vb[1:3, t]) + sum(I3vb[1:3, t])
    
    for (g in 1:3) {
      # Shedding from P compartments (pre-symptomatic)
      shed_P[g, t] <- 
        params$shed_alpha * P[g, t] +
        params$shed_P1a * Pva[g, t] +
        params$shed_P1b * Pvb[g, t]
      
      # Shedding from A compartments (asymptomatic)
      shed_A[g, t] <- 
        params$shed_A1 * A1[g, t] + params$shed_A2 * A2[g, t] + params$shed_A3 * A3[g, t] +
        params$shed_A1a * A1va[g, t] + params$shed_A2a * A2va[g, t] + params$shed_A3a * A3va[g, t] +
        params$shed_A1b * A1vb[g, t] + params$shed_A2b * A2vb[g, t] + params$shed_A3b * A3vb[g, t]
      
      # Shedding from I compartments (symptomatic)
      shed_I[g, t] <- 
        params$shed_I1 * I1[g, t] + params$shed_I2 * I2[g, t] + params$shed_I3 * I3[g, t] +
        params$shed_I1a * I1va[g, t] + params$shed_I2a * I2va[g, t] + params$shed_I3a * I3va[g, t] +
        params$shed_I1b * I1vb[g, t] + params$shed_I2b * I2vb[g, t] + params$shed_I3b * I3vb[g, t]
    }
    
    
    # Total shedding across all groups (optional, if still needed)
    daily_shedding[t] <- sum(shed_P[1:3, t]) + sum(shed_A[1:3, t]) + sum(shed_I[1:3, t])
    
    
    for (g in 1:3) {
      
      # Pre-symptomatic (P)
      P_total[g, t] <- P[g, t] + Pva[g, t] + Pvb[g, t]
      
      # Asymptomatic (A)
      A_total[g, t] <- A1[g, t] + A2[g, t] + A3[g, t] +
        A1va[g, t] + A2va[g, t] + A3va[g, t] +
        A1vb[g, t] + A2vb[g, t] + A3vb[g, t]
      
      # Symptomatic (I)
      I_total[g, t] <- I1[g, t] + I2[g, t] + I3[g, t] +
        I1va[g, t] + I2va[g, t] + I3va[g, t] +
        I1vb[g, t] + I2vb[g, t] + I3vb[g, t]
    }
    
  }
  
  ## Advection-dispersion-delay process
  # Define advection-dispersion-decay using intermediate variables
  
  ## Advection-dispersion-delay process
  sigma <- params$transit_mean * params$transit_cv
  tmax <- ceiling(params$transit_mean + 4 * sigma)
  g <- numeric(tmax)
  for (i in 1:tmax) {
    g[i] <- (1 / (sigma * sqrt(2 * pi))) * exp(-(i - params$transit_mean)^2 / (2 * sigma^2))
  }
  
  contrib <- matrix(0, nrow = tmax, ncol = T)
  delayed_conc <- rep(0, T)
  cp_per_person_mL_all <- rep(0, T)
  log10_conc_all <- rep(0, T)
  
  for (t in 1:T) {
    for (i in 1:tmax) {
      idx <- t - i + 1
      if (idx >= 1) {
        contrib[i, t] <- daily_shedding[idx] * g[i] * exp(-params$mu * i)
      }
    }
    delayed_conc[t] <- sum(contrib[, t])
    cp_per_person_all <- delayed_conc[t] * params$mult / wwtp_population
    cp_per_person_mL_all[t] <- cp_per_person_all * flow_mlalldaily[t]
    log10_conc_all[t] <- log10(cp_per_person_mL_all[t] + 1)  # Avoid log(0)
  }
  
  # ---- Construct symptomatic incidence (I1 transitions) ----
  sympt_new <- colSums(new_I1s + new_I1vas + new_I1vbs)
  
  # ---- Construct cumulative infections (final size proxy) ----
  R_total <- colSums(R)
  
  return(list(
    
    # ---------------------------
    # Surveillance quantities
    # ---------------------------
    cases       = total_new_cases,   # reported cases (model-generated)
    sympt_new   = sympt_new,         # symptomatic incidence (for UF_sympt)
    R_total     = R_total,           # cumulative removed (for UF_total)
    
    # ---------------------------
    # Core epidemic outputs
    # ---------------------------
    prevalence  = active_infected,
    cumulative_incidence = total_Cuminc,
    FOI         = total_lambda,
    
    # ---------------------------
    # Wastewater outputs
    # ---------------------------
    viral_load  = log10_conc_all,
    daily_shedding = daily_shedding,
    
    # ---------------------------
    # Shedding decomposition
    # ---------------------------
    shed_P  = shed_P,
    shed_A  = shed_A,
    shed_I  = shed_I,
    
    # ---------------------------
    # Infection stage totals
    # ---------------------------
    P_total = P_total,
    A_total = A_total,
    I_total = I_total
  ))
}
# ---- Initial conditions based on posterior median ----
# Seed initial infections proportional to contact structure

E_start  <- 14 * norm_contact_dy
P_start  <- 5  * norm_contact_dy
A1_start <- 5  * norm_contact_dy
I1_start <- 5  * norm_contact_dy

A2_start     <- c(0, 0, 0)
A3_start     <- c(0, 0, 0)
I2_start     <- c(0, 0, 0)
I3_start     <- c(0, 0, 0)
R_start      <- c(0, 0, 0)
Cuminc_start <- c(0, 0, 0)

# Susceptible population
S_start <- N - (E_start + P_start + A1_start + I1_start)


# ---- Initial conditions: first-dose vaccinated ----
Va_start   <- c(0, 0, 0)
Eva_start  <- c(0, 0, 0)
Pva_start  <- c(0, 0, 0)
A1va_start <- c(0, 0, 0)
A2va_start <- c(0, 0, 0)
A3va_start <- c(0, 0, 0)
I1va_start <- c(0, 0, 0)
I2va_start <- c(0, 0, 0)
I3va_start <- c(0, 0, 0)


# ---- Initial conditions: second-dose vaccinated ----
Vb_start   <- c(0, 0, 0)
Evb_start  <- c(0, 0, 0)
Pvb_start  <- c(0, 0, 0)
A1vb_start <- c(0, 0, 0)
A2vb_start <- c(0, 0, 0)
A3vb_start <- c(0, 0, 0)
I1vb_start <- c(0, 0, 0)
I2vb_start <- c(0, 0, 0)
I3vb_start <- c(0, 0, 0)


# ---- Vaccination coverage inputs ----
case_dat2 <- read_excel("data/raw/case_data_V2.xlsx", sheet = "Mergeddat")

prop_cov  <- case_dat2$cov_dose1
prop_cov2 <- case_dat2$cov_dose2


# ---- Convert first-dose coverage to daily per-capita rate ----
p <- prop_cov
v <- rep(NA, length(p))

v[1] <- -log(1 - p[1])

for (i in 2:length(p)) {
  x <- (1 - sum(p[1:i])) * exp(sum(v[1:(i - 1)]))
  v[i] <- -log(x)
}

# prepend burn-in period
alpha_v <- c(rep(0, burn_in_timesteps), v)


# ---- Convert second-dose coverage to daily per-capita rate ----
p2 <- prop_cov2
v2 <- rep(NA, length(p2))

v2[1] <- -log(1 - p2[1])

for (i in 2:length(p2)) {
  x2 <- (1 - sum(p2[1:i])) * exp(sum(v2[1:(i - 1)]))
  v2[i] <- -log(x2)
}

# prepend burn-in period
mu_v <- c(rep(0, burn_in_timesteps), v2)


# ---- Assemble initial state list ----
init_vals <- list(
  # unvaccinated compartments
  S = S_start,
  E = E_start,
  P = P_start,
  A1 = A1_start,
  A2 = A2_start,
  A3 = A3_start,
  I1 = I1_start,
  I2 = I2_start,
  I3 = I3_start,
  R = R_start,
  Cuminc = Cuminc_start,
  
  # first-dose vaccinated compartments
  Va = Va_start,
  Eva = Eva_start,
  Pva = Pva_start,
  A1va = A1va_start,
  A2va = A2va_start,
  A3va = A3va_start,
  I1va = I1va_start,
  I2va = I2va_start,
  I3va = I3va_start,
  
  # second-dose vaccinated compartments
  Vb = Vb_start,
  Evb = Evb_start,
  Pvb = Pvb_start,
  A1vb = A1vb_start,
  A2vb = A2vb_start,
  A3vb = A3vb_start,
  I1vb = I1vb_start,
  I2vb = I2vb_start,
  I3vb = I3vb_start,
  
  # vaccination rates
  alpha_v = alpha_v,
  mu_v = mu_v
)


# ---- Reduced shedding scenarios ----
remaining_fracs <- c(0.8, 0.6, 0.2)

runs <- map(remaining_fracs, function(f) {
  
  out <- mysimulate_model(
    params, init_vals, contact, N, T,
    vax_scenario = "reduced_shedding",
    v1_frac = f,
    v2_frac = f
  )
  
  out$label <- paste0(
    "Reduced shedding: ",
    round((1 - f) * 100),
    "%"
  )
  
  out
})
# ==================================================
# Scenario simulations
# ==================================================

# ---- 1. Baseline scenario ----
baseline_result <- list(
  "Baseline" = mysimulate_model(
    params, init_vals, contact, N, T,
    vax_scenario = "baseline"
  )
)


# ---- 2. Reduced shedding scenarios ----
remaining_fracs <- c(0.8, 0.6, 0.2)

runs <- map(remaining_fracs, function(f) {
  
  out <- mysimulate_model(
    params, init_vals, contact, N, T,
    vax_scenario = "reduced_shedding",
    v1_frac = f,
    v2_frac = f
  )
  
  out$label <- paste0(
    "Reduced shedding: ",
    round((1 - f) * 100),
    "%"
  )
  
  out
})


# ---- 3. Transit time scenarios (posterior-anchored) ----
transit_results <- run_transit_scenarios(
  params, init_vals, contact, N, T,
  factors = transit_factors
)


# ---- 4. Asymptomatic proportion scenarios ----
asym_results <- run_asym_scenarios(
  params, init_vals, contact, N, T,
  m_factors
)


# ---- 5. Reporting delay scenarios ----
reporting_delays <- c(5, 10, 14)

reporting_delay_results <- run_reporting_delay_scenarios(
  params, init_vals, contact, N, T,
  delays = reporting_delays
)


# ---- 6. Combine all scenarios into master list ----
scenario_list <- c(
  baseline_result,
  setNames(
    runs,
    paste0(
      "Reduced shedding: ",
      round((1 - remaining_fracs) * 100),
      "%"
    )
  ),
  transit_results,
  asym_results,
  reporting_delay_results
)
# ==================================================
# Posterior sampling for counterfactual simulations
# ==================================================

# ---- Number of posterior draws ----
n_draws <- 200        # increase (e.g. 5000) for full analysis
G <- 3                # number of activity groups
T <- 199              # number of time points

draw_indices <- sample(
  1:nrow(mcmc_matrixallcom),
  size = n_draws,
  replace = FALSE
)


# ---- Reduced shedding fractions ----
fractions <- c(0.8, 0.6, 0.2)


# ---- Storage objects ----
# Structure:
#   - baseline
#   - reduced_shedding (one entry per fraction)

results <- list(
  
  baseline = list(
    cases  = matrix(NA, n_draws, T),
    vload  = matrix(NA, n_draws, T),
    cuminc = matrix(NA, n_draws, T),
    prev   = matrix(NA, n_draws, T),
    FOI    = matrix(NA, n_draws, T),
    sympt_new = matrix(NA, n_draws, T),
    R_total   = matrix(NA, n_draws, T),
    
    shed_P  = array(NA, dim = c(n_draws, G, T)),
    shed_A  = array(NA, dim = c(n_draws, G, T)),
    shed_I  = array(NA, dim = c(n_draws, G, T)),
    
    P_total = array(NA, dim = c(n_draws, G, T)),
    A_total = array(NA, dim = c(n_draws, G, T)),
    I_total = array(NA, dim = c(n_draws, G, T))
  ),
  
  reduced_shedding = vector("list", length(fractions))
)


# ---- Pre-allocate storage for each reduced shedding scenario ----
for (f_idx in seq_along(fractions)) {
  
  results$reduced_shedding[[f_idx]] <- list(
    frac   = fractions[f_idx],
    
    cases  = matrix(NA, n_draws, T),
    vload  = matrix(NA, n_draws, T),
    cuminc = matrix(NA, n_draws, T),
    prev   = matrix(NA, n_draws, T),
    FOI    = matrix(NA, n_draws, T),
    sympt_new = matrix(NA, n_draws, T),
    R_total   = matrix(NA, n_draws, T),
    
    shed_P  = array(NA, dim = c(n_draws, G, T)),
    shed_A  = array(NA, dim = c(n_draws, G, T)),
    shed_I  = array(NA, dim = c(n_draws, G, T)),
    
    P_total = array(NA, dim = c(n_draws, G, T)),
    A_total = array(NA, dim = c(n_draws, G, T)),
    I_total = array(NA, dim = c(n_draws, G, T))
  )
}

# Extend results container: add a list for transit scenarios
transit_factors <-transit_factors

# Extend results container: add a list for transit scenarios
results$longer_transit <- vector("list", length(transit_factors))

for (tf_idx in seq_along(transit_factors)) {
  results$longer_transit[[tf_idx]] <- list(
    factor = transit_factors[tf_idx],
    cases = matrix(NA, n_draws, T),
    vload = matrix(NA, n_draws, T),
    cuminc = matrix(NA, n_draws, T),
    prev = matrix(NA, n_draws, T),
    FOI = matrix(NA, n_draws, T),
    sympt_new = matrix(NA, n_draws, T),
    R_total   = matrix(NA, n_draws, T),
    shed_P = array(NA, dim = c(n_draws, G, T)),
    shed_A = array(NA, dim = c(n_draws, G, T)),
    shed_I = array(NA, dim = c(n_draws, G, T)),
    P_total = array(NA, dim = c(n_draws, G, T)),
    A_total = array(NA, dim = c(n_draws, G, T)),
    I_total = array(NA, dim = c(n_draws, G, T))
  )
}


# Asymptomatic proportion multipliers
m_factors <- c(1.20, 1.60, 1.80)

# Extend results container for asymptomatic scenarios
results$asym_prop <- vector("list", length(m_factors))

for (mf_idx in seq_along(m_factors)) {
  results$asym_prop[[mf_idx]] <- list(
    factor = m_factors[mf_idx],
    cases = matrix(NA, n_draws, T),
    vload = matrix(NA, n_draws, T),
    cuminc = matrix(NA, n_draws, T),
    prev = matrix(NA, n_draws, T),
    FOI = matrix(NA, n_draws, T),
    sympt_new = matrix(NA, n_draws, T),
    R_total   = matrix(NA, n_draws, T),
    shed_P = array(NA, dim = c(n_draws, G, T)),
    shed_A = array(NA, dim = c(n_draws, G, T)),
    shed_I = array(NA, dim = c(n_draws, G, T)),
    P_total = array(NA, dim = c(n_draws, G, T)),
    A_total = array(NA, dim = c(n_draws, G, T)),
    I_total = array(NA, dim = c(n_draws, G, T))
  )
}

# Reporting delay scenarios

results$reporting_delay <- vector("list", length(reporting_delays))
for (rd_idx in seq_along(reporting_delays)) {
  results$reporting_delay[[rd_idx]] <- list(
    delay = reporting_delays[rd_idx],
    cases = matrix(NA, n_draws, T),
    vload = matrix(NA, n_draws, T),
    cuminc = matrix(NA, n_draws, T),
    prev = matrix(NA, n_draws, T),
    FOI = matrix(NA, n_draws, T),
    sympt_new = matrix(NA, n_draws, T),
    R_total   = matrix(NA, n_draws, T),
    shed_P = array(NA, dim = c(n_draws, G, T)),
    shed_A = array(NA, dim = c(n_draws, G, T)),
    shed_I = array(NA, dim = c(n_draws, G, T)),
    P_total = array(NA, dim = c(n_draws, G, T)),
    A_total = array(NA, dim = c(n_draws, G, T)),
    I_total = array(NA, dim = c(n_draws, G, T))
  )
}


# Loop over posterior draws
for (i in seq_along(draw_indices)) {
  draw <- draw_indices[i]
  params <- extract_params(mcmc_matrixallcom[draw, ])
  
  # ---- Baseline ----
  sim_base <- mysimulate_model(params, init_vals, contact, N, T,
                               vax_scenario = "baseline")
  
  results$baseline$cases[i, ]    <- sim_base$cases
  results$baseline$vload[i, ]    <- sim_base$viral_load
  results$baseline$cuminc[i, ]   <- sim_base$cumulative_incidence
  results$baseline$prev[i, ]     <- sim_base$prevalence
  results$baseline$FOI[i, ]      <- sim_base$FOI
  # ---- ADD THESE TWO LINES ----
  results$baseline$sympt_new[i, ] <- sim_base$sympt_new
  results$baseline$R_total[i, ]   <- sim_base$R_total
  results$baseline$shed_P[i, , ] <- sim_base$shed_P
  results$baseline$shed_A[i, , ] <- sim_base$shed_A
  results$baseline$shed_I[i, , ] <- sim_base$shed_I
  results$baseline$P_total[i, , ] <- sim_base$P_total
  results$baseline$A_total[i, , ] <- sim_base$A_total
  results$baseline$I_total[i, , ] <- sim_base$I_total
  
  # ---- Reduced shedding scenarios ----
  for (f_idx in seq_along(fractions)) {
    f <- fractions[f_idx]
    sim_rs <- mysimulate_model(params, init_vals, contact, N, T,
                               vax_scenario = "reduced_shedding",
                               v1_frac = f, v2_frac = f)
    
    rs <- results$reduced_shedding[[f_idx]]
    rs$cases[i, ]     <- sim_rs$cases
    rs$vload[i, ]     <- sim_rs$viral_load
    rs$cuminc[i, ]    <- sim_rs$cumulative_incidence
    rs$prev[i, ]      <- sim_rs$prevalence
    rs$FOI[i, ]       <- sim_rs$FOI
    # ---- ADD THESE TWO LINES ----
    rs$sympt_new[i, ] <- sim_rs$sympt_new
    rs$R_total[i, ]   <- sim_rs$R_total
    rs$shed_P[i, , ]  <- sim_rs$shed_P
    rs$shed_A[i, , ]  <- sim_rs$shed_A
    rs$shed_I[i, , ]  <- sim_rs$shed_I
    rs$P_total[i, , ] <- sim_rs$P_total
    rs$A_total[i, , ] <- sim_rs$A_total
    rs$I_total[i, , ] <- sim_rs$I_total
    
    # Save back into results list
    results$reduced_shedding[[f_idx]] <- rs
  }
  
  # ---- NEW: Longer transit scenarios (change only transit_mean) ----
  for (tf_idx in seq_along(transit_factors)) {
    tf <- transit_factors[tf_idx]
    p2 <- params
    p2$transit_mean <- params$transit_mean * tf  # keep transit_cv as is
    
    sim_tr <- mysimulate_model(p2, init_vals, contact, N, T, vax_scenario = "baseline")
    
    tr <- results$longer_transit[[tf_idx]]
    tr$cases[i, ]     <- sim_tr$cases          # same epi, but we save for completeness
    tr$vload[i, ]     <- sim_tr$viral_load     # WW will shift/smooth
    tr$cuminc[i, ]    <- sim_tr$cumulative_incidence
    tr$prev[i, ]      <- sim_tr$prevalence
    tr$FOI[i, ]       <- sim_tr$FOI
    # ---- ADD THESE TWO LINES ----
    tr$sympt_new[i, ] <- sim_tr$sympt_new
    tr$R_total[i, ]   <- sim_tr$R_total
    tr$shed_P[i, , ]  <- sim_tr$shed_P
    tr$shed_A[i, , ]  <- sim_tr$shed_A
    tr$shed_I[i, , ]  <- sim_tr$shed_I
    tr$P_total[i, , ] <- sim_tr$P_total
    tr$A_total[i, , ] <- sim_tr$A_total
    tr$I_total[i, , ] <- sim_tr$I_total
    results$longer_transit[[tf_idx]] <- tr
  }
  
  
  # ---- NEW: Asymptomatic proportion scenarios (scale m) ----
  for (mf_idx in seq_along(m_factors)) {
    mf <- m_factors[mf_idx]
    p3 <- params
    p3$m <- p3$m * mf  # scale posterior-estimated m
    
    sim_m <- mysimulate_model(p3, init_vals, contact, N, T, vax_scenario = "baseline")
    
    am <- results$asym_prop[[mf_idx]]
    am$cases[i, ]     <- sim_m$cases
    am$vload[i, ]     <- sim_m$viral_load
    am$cuminc[i, ]    <- sim_m$cumulative_incidence
    am$prev[i, ]      <- sim_m$prevalence
    am$FOI[i, ]       <- sim_m$FOI
    # ---- ADD THESE TWO LINES ----
    am$sympt_new[i, ] <- sim_m$sympt_new
    am$R_total[i, ]   <- sim_m$R_total
    
    am$shed_P[i, , ]  <- sim_m$shed_P
    am$shed_A[i, , ]  <- sim_m$shed_A
    am$shed_I[i, , ]  <- sim_m$shed_I
    am$P_total[i, , ] <- sim_m$P_total
    am$A_total[i, , ] <- sim_m$A_total
    am$I_total[i, , ] <- sim_m$I_total
    
    results$asym_prop[[mf_idx]] <- am
  }
  
  # ---- NEW: Reporting delay scenarios ----
  # ---- NEW: Reporting delay scenarios ----
  for (rd_idx in seq_along(reporting_delays)) {
    rd <- reporting_delays[rd_idx]
    
    p4 <- params
    p4$report_lag <- rd  # override default 7-day lag
    
    sim_rd <- mysimulate_model(
      p4, init_vals, contact, N, T,
      vax_scenario = "baseline"
    )
    
    rd_res <- results$reporting_delay[[rd_idx]]
    rd_res$cases[i, ]     <- sim_rd$cases
    rd_res$vload[i, ]     <- sim_rd$viral_load
    rd_res$cuminc[i, ]    <- sim_rd$cumulative_incidence
    rd_res$prev[i, ]      <- sim_rd$prevalence
    rd_res$FOI[i, ]       <- sim_rd$FOI
    # ---- ADD THESE TWO LINES ----
    rd_res$sympt_new[i, ] <- sim_rd$sympt_new
    rd_res$R_total[i, ]   <- sim_rd$R_total
    rd_res$shed_P[i, , ]  <- sim_rd$shed_P
    rd_res$shed_A[i, , ]  <- sim_rd$shed_A
    rd_res$shed_I[i, , ]  <- sim_rd$shed_I
    rd_res$P_total[i, , ] <- sim_rd$P_total
    rd_res$A_total[i, , ] <- sim_rd$A_total
    rd_res$I_total[i, , ] <- sim_rd$I_total
    
    results$reporting_delay[[rd_idx]] <- rd_res
  }
}

# ==================================================
# Under-ascertainment calculations
# ==================================================

# Observed reported cases
obs_cases <- sum(cases_obsb)
burn <- burn_in_timesteps


# ----------------------------
# Helper function
# ----------------------------
compute_UF <- function(sympt_mat, R_mat, burn, obs_cases) {
  
  total_sympt <- rowSums(sympt_mat[, -(1:burn), drop = FALSE])
  UF_sympt <- total_sympt / obs_cases
  
  final_size <- R_mat[, ncol(R_mat)]
  UF_total <- final_size / obs_cases
  
  list(
    UF_sympt = c(
      median   = median(UF_sympt),
      lower_95 = quantile(UF_sympt, 0.025),
      upper_95 = quantile(UF_sympt, 0.975)
    ),
    UF_total = c(
      median   = median(UF_total),
      lower_95 = quantile(UF_total, 0.025),
      upper_95 = quantile(UF_total, 0.975)
    )
  )
}


# ----------------------------
# Baseline
# ----------------------------
UF_baseline <- compute_UF(
  results$baseline$sympt_new,
  results$baseline$R_total,
  burn,
  obs_cases
)


# ----------------------------
# Reduced shedding scenarios
# ----------------------------
UF_reduced <- lapply(results$reduced_shedding, function(rs) {
  compute_UF(rs$sympt_new, rs$R_total, burn, obs_cases)
})


# ----------------------------
# Transit scenarios
# ----------------------------
UF_transit <- lapply(results$longer_transit, function(tr) {
  compute_UF(tr$sympt_new, tr$R_total, burn, obs_cases)
})


# ----------------------------
# Asymptomatic proportion scenarios
# ----------------------------
UF_asym <- lapply(results$asym_prop, function(am) {
  compute_UF(am$sympt_new, am$R_total, burn, obs_cases)
})


# ----------------------------
# Reporting delay scenarios
# ----------------------------
UF_reporting <- lapply(results$reporting_delay, function(rd) {
  compute_UF(rd$sympt_new, rd$R_total, burn, obs_cases)
})


sapply(results$reduced_shedding, function(rs) {
  median(rowSums(rs$sympt_new[, -(1:burn), drop = FALSE]))
})

sapply(results$reduced_shedding, function(rs) {
  median(rs$R_total[, ncol(rs$R_total)])
})




