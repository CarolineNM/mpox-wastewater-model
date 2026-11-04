set.seed(123)  # reproducible MCMC initialization

# ==================================================
# Script: 03_fit_combined.R
# Purpose: Fit combined case + wastewater mpox transmission model
# Project: Mpox_clean
# Updated: 2026-04-22
# ==================================================

# ---- Packages ----
library(R2jags)
library(runjags)
library(mcmcplots)
library(tidyverse)
library(readxl)
library(EnvStats)

options(scipen = 999)

# ---- Paths ----
data_file <- file.path("data", "raw", "case_data_v2.xlsx")
out_dir   <- file.path("outputs", "fits")

# ---- Load data ----
case_dat <- read_excel(data_file, sheet = "cases")
cov_dat  <- read_excel(data_file, sheet = "Mergeddat")
ww_dat   <- read_excel(data_file, sheet = "dailyWW")
ww_flow  <- read_excel(data_file, sheet=  "flowrate")

# ---- Model setup ----

textstring <- "
model {

  # ---- Wastewater observation parameters ----
  mu <- mu_fixed

  log_mult ~ dnorm(log(1.5e-4), 0.05) T(log(1e-8), log(1e0))
  mult <- exp(log_mult)

  tau_ww ~ dgamma(40, 48)

  transit_time_mean ~ dnorm(2.5, 0.25) T(0.1, 10)
  transit_time_cv   ~ dnorm(0.3, 3) T(0.1, 1)

  # ---- Transmission and mixing ----
  beta ~ dnorm(0.8, 100) T(0, 1)

  # assortative mixing prior; mean ~0.95
  kappa ~ dbeta(40, 2)

  # negative binomial overdispersion prior; mean = 4
  phi ~ dgamma(2, 0.5)


  # ---- Initial seeding conditions ----
  E0  ~ dpois(5)
  P0  ~ dpois(1)
  A10 ~ dpois(1)
  I10 ~ dpois(1)


  # ---- Exposed period (E compartment) ----
  # mean duration ~6.26 days
  delta_inv ~ dgamma(300.65298, 47.98256)

  delta   <- 1 / delta_inv
  deltava <- delta
  deltavb <- delta


  # ---- Pre-symptomatic period (P compartment) ----
  # fixed duration = 2 days
  epsilon_inv <- 2

  epsilon   <- 1 / epsilon_inv
  epsilonva <- epsilon
  epsilonvb <- epsilon


  # ---- Probability asymptomatic ----
  # mean ~34%
  m ~ dbeta(49.09, 94.88)

  mva <- m
  mvb <- m


  # ---- Symptomatic infectious period (I1-I3) ----
  # total duration ~21 days across three compartments
  theta_invall ~ dgamma(42.3, 1.99)

  theta_inv <- theta_invall / 3
  theta     <- 1 / theta_inv

  thetava <- theta
  thetavb <- theta


  # ---- Asymptomatic infectious period (A1-A3) ----
  # total duration ~14 days across three compartments
  omega_invall ~ dgamma(20.6, 1.45)

  omega_inv <- omega_invall / 3
  omega     <- 1 / omega_inv

  omegava <- omega
  omegavb <- omega


  # ---- Reporting process ----
  # prior mean = 0.50
  report_frac ~ dbeta(10, 10)

  # fixed reporting delay (days)
  report_lag <- 7


  # ---- Vaccine effectiveness: first dose ----
  # mean effectiveness ~35.8%
  Vea ~ dbeta(49.3, 87.4)

  # residual susceptibility after dose 1
  rsa <- 1 - Vea


  # ---- Vaccine effectiveness: second dose ----
  # mean effectiveness ~66%
  Veb ~ dbeta(69.3, 35.6)

  # residual susceptibility after dose 2
  rsb <- 1 - Veb

  # Initial conditions for each group
 for (g in 1:3) {
    E_start[g] <- E0 * norm_contact_dy[g]
    P_start[g] <- P0 * norm_contact_dy[g]
    A1_start[g] <- A10 * norm_contact_dy[g]
    I1_start[g] <- I10 * norm_contact_dy[g]
    E[g,1]<-E_start[g]
    P[g, 1] <- P_start[g]
    A1[g, 1] <- A1_start[g]
    I1[g, 1] <- I1_start[g]
    S[g, 1] <- max(0,N[g] - E_start[g]-P_start[g]-A1_start[g]-I1_start[g])
    A2[g, 1] <- A2_start[g]
    A3[g, 1] <- A3_start[g]
    I2[g, 1] <- I2_start[g]
    I3[g, 1] <- I3_start[g]
    R[g, 1] <- R_start[g]
    Cuminc[g, 1]<-Cuminc_start[g]
    
    #Inits for vaccinated comps va
    Va[g,1]<-Va_start[g]
    Eva[g,1]<-Eva_start[g]
    Pva[g, 1] <- Pva_start[g]
    A1va[g, 1] <- A1va_start[g]
    A2va[g, 1] <- A2va_start[g]
    A3va[g, 1] <- A3va_start[g]
    I1va[g, 1] <- I1va_start[g]
    I2va[g, 1] <- I2va_start[g]
    I3va[g, 1] <- I3va_start[g]
    
    #Inits for vaccinated comps vb
    Vb[g,1]<-Vb_start[g]
    Evb[g,1]<-Evb_start[g]
    Pvb[g, 1] <- Pvb_start[g]
    A1vb[g, 1] <- A1vb_start[g]
    A2vb[g, 1] <- A2vb_start[g]
    A3vb[g, 1] <- A3vb_start[g]
    I1vb[g, 1] <- I1vb_start[g]
    I2vb[g, 1] <- I2vb_start[g]
    I3vb[g, 1] <- I3vb_start[g]
    
    for (d in 1:va_delay) {
    va_queue[g,d,1] <- 0
    vb_queue[g,d,1] <- 0
  }
 }

  # Define the assortative mixing matrix once outside of the loop
  for (g in 1:3) {
    for (g_ in 1:3) {
      mix_ass[g, g_] <- ifelse(g == g_, 1, 0)
    }
  }
  
  # Define the proportional mixing matrix once outside of the loop
  for (g in 1:3) {
    for (g_ in 1:3) {
      mix_prp[g, g_] <- (contact[g_] * N[g_]) / sum(contact * N)
    }
  }
  
  # Initialize prevalence for t = 1
  for (g in 1:3) {
    prevalence[g, 1] <- ((tau*P[g, 1]) + I1[g, 1] + I2[g, 1] + I3[g, 1]+
                          (tau*A1[g, 1]) + (tau*A2[g, 1]) + (tau*A3[g, 1])+
                          (tau*Pva[g, 1])+I1va[g, 1] + I2va[g, 1] + I3va[g, 1]+
                          (tau*A1va[g, 1])+ (tau*A2va[g, 1])+ (tau*A3va[g, 1])+
                          (tau*Pvb[g, 1])+I1vb[g, 1] + I2vb[g, 1] + I3vb[g, 1]+
                          (tau*A1vb[g, 1])+ (tau*A2vb[g, 1])+ (tau*A3vb[g, 1]))/ N[g]
  }
  

  # Initialize mixing matrix (mix) for t = 1
  for (g in 1:3) {
    for (g_ in 1:3) {
      mix[g, g_, 1] <- kappa * mix_ass[g, g_] + (1 - kappa) * mix_prp[g, g_]
    }
  }

 
 # Calculate initial lambda_det for t = 1 based on the initial conditions
for (g in 1:3) {
  lambda_det[g, 1] <- sum(beta * contact[g] * mix[g, 1:3, 1] * prevalence[1:3, 1])
}

total_lambda[1] <- sum(lambda_det[1:3, 1])

 #calculating new cases to include in likelihood

    total_new_cases[1] <- 0                    # Initial people in Icomp
    total_Cuminc[1] <- sum(Cuminc[1:3, 1])     # Initial cumulative Incidence
    
     active_infected[1] <- 0                  # Initial prevalence
     daily_shedding[1]<-0                     #Initial shed viral load
      infectious_pop[1]<-0
     shed_per_inf[1]<-0
     sympt_new[1]<-0
     R_total[1]<-0
     
     # Initial group-specific daily shedding at t = 1
for (g in 1:3) {
  daily_shedding_g[g, 1] <- 0
}


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
          new_reported_cases[g, t] <- step(t - report_lag) * report_frac *
                            (new_I1s[g, max(1, t - report_lag)]+new_I1vas[g, max(1, t - report_lag)]+
                            new_I1vbs[g, max(1, t - report_lag)])

    }
    
       # Total new reported cases at time t(incidence at t)
        total_new_cases[t] <- sum(new_reported_cases[1:3, t])
        
        total_Cuminc[t] <- sum(Cuminc[1:3, t])
    
        
         # total symptomatic incidence
    sympt_new[t] <- new_I1s[1, t-1] + new_I1s[2, t-1] + new_I1s[3, t-1] +
                    new_I1vas[1, t-1] + new_I1vas[2, t-1] + new_I1vas[3, t-1] +
                    new_I1vbs[1, t-1] + new_I1vbs[2, t-1] + new_I1vbs[3, t-1]

    # total removed (cumulative infections)
    R_total[t] <- R[1, t] + R[2, t] + R[3, t]
        
        
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
      
      #Shedding Parameters (Fixed log10 converted outside)
      ##Viral Shedding Calculation (Daily)
      
      for (g in 1:3) {

  # 1. Shedding from P compartments (pre-symptomatic)
  shed_P[g, t] <- shed_alpha * (
    P[g, t] +
    Pva[g, t] +
    Pvb[g, t]
  )

  # 2. Shedding from A compartments (asymptomatic)
  shed_A[g, t] <- 
    shed_A1 * A1[g, t] + shed_A2 * A2[g, t] + shed_A3 * A3[g, t] +
    shed_A1a * A1va[g, t] + shed_A2a * A2va[g, t] + shed_A3a * A3va[g, t] +
    shed_A1b * A1vb[g, t] + shed_A2b * A2vb[g, t] + shed_A3b * A3vb[g, t]

  # 3. Shedding from I compartments (symptomatic)
  shed_I[g, t] <- 
    shed_I1 * I1[g, t] + shed_I2 * I2[g, t] + shed_I3 * I3[g, t] +
    shed_I1a * I1va[g, t] + shed_I2a * I2va[g, t] + shed_I3a * I3va[g, t] +
    shed_I1b * I1vb[g, t] + shed_I2b * I2vb[g, t] + shed_I3b * I3vb[g, t]
}

# Total shedding across all groups (optional, if still needed)
daily_shedding[t] <- sum(shed_P[1:3, t]) + sum(shed_A[1:3, t]) + sum(shed_I[1:3, t])

# Group-specific daily shedding
for (g in 1:3) {
  daily_shedding_g[g, t] <- shed_P[g, t] + shed_A[g, t] + shed_I[g, t]
}


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
  
  infectious_pop[t] <- sum(P_total[1:3,t]) + sum(A_total[1:3,t]) + sum(I_total[1:3,t])
  shed_per_inf[t] <- daily_shedding[t] / infectious_pop[t]

}
  
  ## Advection-dispersion-delay process
  # Define advection-dispersion-decay using intermediate variables

sigma <- transit_time_mean * transit_time_cv  

for (i in 1:tmax) {
  g[i] <- (1 / (sigma * sqrt(2 * 3.14159))) * exp(-pow(i - transit_time_mean, 2) / (2 * pow(sigma, 2)))
}

for (t in 1:T) {
  for (i in 1:tmax) {
    lag[i, t] <- t - i + 1
    is_valid[i, t] <- step(lag[i, t])
    safe_index[i, t] <- max(lag[i, t], 1)
    safe_shedding[i, t] <- is_valid[i, t] * daily_shedding[safe_index[i, t]]
    contrib[i, t] <- safe_shedding[i, t] * g[i] * exp(-mu * i)
  }
  delayed_conc[t] <- sum(contrib[1:tmax, t])
}

for (i in 1:tmax) {
  g_kernel[i] <- (1 / (sigma * sqrt(2 * 3.14159))) * exp(-pow(i - transit_time_mean, 2) / (2 * pow(sigma, 2)))
}

# Group-specific delayed concentration
for (g in 1:3) {
  for (t in 1:T) {
    for (i in 1:tmax) {
      lag_g[i, t, g] <- t - i + 1
      is_valid_g[i, t, g] <- step(lag_g[i, t, g])
      safe_index_g[i, t, g] <- max(lag_g[i, t, g], 1)
      safe_shedding_g[i, t, g] <- is_valid_g[i, t, g] * daily_shedding_g[g, safe_index_g[i, t, g]]
      contrib_g[i, t, g] <- safe_shedding_g[i, t, g] * g_kernel[i] * exp(-mu * i)
    }
    delayed_conc_g[g, t] <- sum(contrib_g[1:tmax, t, g])
  }
}


  ##Case likelihood
for (t in (burn_in_timesteps + 1):(T_caseobs + burn_in_timesteps)) {
    mu_nb[t] <- total_new_cases[t]+ 1e-6 # Mean for the Negative Binomial to ensure its not neg
    p_nb[t]<-phi/(phi+mu_nb[t]) ##convert dispersion param to success prob
    cases_obs[t - burn_in_timesteps] ~ dnegbin(p_nb[t], phi)  # Negative Binomial likelihood
      # Posterior predictive distribution
    cases_pred[t - burn_in_timesteps] ~ dnegbin(p_nb[t], phi)
}

 ##Viral load likelihood
  ##Scaling, normalization and log transformation
  
 for (t in (burn_in_timesteps + 1):(T_caseobs + burn_in_timesteps)) {
  cp_total_all[t] <- delayed_conc[t] * mult
  cp_per_person_all[t] <- cp_total_all[t] / wwtp_population
  cp_per_person_mL_all[t] <- cp_per_person_all[t] * flow_mlalldaily[t - burn_in_timesteps]
  log10_conc_all[t] <- log(cp_per_person_mL_all[t] + 1) / log(10)
}
for (w in 1:T_WWobs) {
   cp_total[w] <- delayed_conc[ww_sample_days[w]] * mult
  cp_per_person[w] <- cp_total[w] / wwtp_population
  #cp_per_person_mL[w] <- cp_per_person[w] * flow_mL_daily[w]
  #log10_conc[w] <- log(cp_per_person_mL[w] + 1) / log(10)
  cp_raw[w] <- cp_per_person[w] * flow_mL_daily[w]
  cp_safe[w] <- max(cp_raw[w], 1e-6)
  log10_conc[w] <- log(cp_safe[w]) / log(10)
  ww_obs[w] ~ dnorm(log10_conc[w], tau_ww)
  ww_pred[w] ~ dnorm(log10_conc[w], tau_ww)

}

# All days (for plotting/PPD over full timeline)
for (g in 1:3) {
  for (t in (burn_in_timesteps + 1):(T_caseobs + burn_in_timesteps)) {
    cp_total_all_g[g, t] <- delayed_conc_g[g, t] * mult
    cp_per_person_all_g[g, t] <- cp_total_all_g[g, t] / wwtp_population
    cp_per_person_mL_all_g[g, t] <- cp_per_person_all_g[g, t] * flow_mlalldaily[t - burn_in_timesteps]
    log10_conc_all_g[g, t] <- log(cp_per_person_mL_all_g[g, t] + 1) / log(10)
  }
}

for (g in 1:3) {
  for (t in (burn_in_timesteps + 1):(T_caseobs + burn_in_timesteps)) {
    # Posterior predictive draw for full daily series (all timesteps)
    log10_conc_all_g_pred[g, t] ~ dnorm(log10_conc_all_g[g, t], tau_ww)
  }
}

for (g in 1:3) {
  for (t in (burn_in_timesteps + 1):(T_caseobs + burn_in_timesteps)) {
    mu_nb_g[g, t] <- new_reported_cases[g, t] + 1e-6
    p_nb_g[g, t] <- phi / (phi + mu_nb_g[g, t])
    cases_pred_g[g, t - burn_in_timesteps] ~ dnegbin(p_nb_g[g, t], phi)
  }}
}"


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

# normalized daily contact rates for distributing initial exposure risk
norm_contact_dy <- contact_dy / sum(contact_dy)


# ---- Population size and simulation horizon ----
N <- c(15660, 9135, 1305)

burn_in_timesteps <- 30
T <- 169 + burn_in_timesteps
T_caseobs <- length(case_dat$total_cases)


# ---- Initial conditions: unvaccinated compartments ----
A2_start <- c(0, 0, 0)
A3_start <- c(0, 0, 0)
I2_start <- c(0, 0, 0)
I3_start <- c(0, 0, 0)
R_start <- c(0, 0, 0)
Cuminc_start <- c(0, 0, 0)


# ---- Initial conditions: first-dose vaccinated compartments ----
Va_start   <- c(0, 0, 0)
Pva_start  <- c(0, 0, 0)
Eva_start  <- c(0, 0, 0)
A1va_start <- c(0, 0, 0)
A2va_start <- c(0, 0, 0)
A3va_start <- c(0, 0, 0)
I1va_start <- c(0, 0, 0)
I2va_start <- c(0, 0, 0)
I3va_start <- c(0, 0, 0)


# ---- Initial conditions: second-dose vaccinated compartments ----
Vb_start   <- c(0, 0, 0)
Pvb_start  <- c(0, 0, 0)
Evb_start  <- c(0, 0, 0)
A1vb_start <- c(0, 0, 0)
A2vb_start <- c(0, 0, 0)
A3vb_start <- c(0, 0, 0)
I1vb_start <- c(0, 0, 0)
I2vb_start <- c(0, 0, 0)
I3vb_start <- c(0, 0, 0)


# ---- Observed case data ----
cases_obsb <- case_dat$total_cases


# ---- Vaccination coverage inputs ----
prop_cov  <- cov_dat$cov_dose1
prop_cov2 <- cov_dat$cov_dose2


# ---- Convert first-dose coverage to daily per-capita vaccination rate ----
p <- prop_cov
v <- rep(NA, length(p))

v[1] <- -log(1 - p[1])

for (i in 2:length(p)) {
  x <- (1 - sum(p[1:i])) * exp(sum(v[1:(i - 1)]))
  v[i] <- -log(x)
}

# prepend burn-in period with zero vaccination
alpha_v <- c(rep(0, burn_in_timesteps), v)


# ---- Convert second-dose coverage to daily per-capita vaccination rate ----
p2 <- prop_cov2
v2 <- rep(NA, length(p2))

v2[1] <- -log(1 - p2[1])

for (i in 2:length(p2)) {
  x2 <- (1 - sum(p2[1:i])) * exp(sum(v2[1:(i - 1)]))
  v2[i] <- -log(x2)
}

# prepend burn-in period with zero vaccination
mu_v <- c(rep(0, burn_in_timesteps), v2)


# ---- Shedding parameters ----
tau <- 0.776

# Stage-specific shedding, log10 copies per mL of viral transport media,
# read from Suner et al. 2023 figure 3A (median viral load by time from
# symptom onset) for lesion and rectal swabs at the bins 1-7, 8-14 and 15-21
# days, which correspond to stages I1, I2 and I3. The two sites are summed on
# the LINEAR scale. Derivation in scripts/derive_shedding_rates.R.
#
# These are transport-media concentrations and cannot be converted to copies
# per person per day, so the absolute scale is not identifiable here and is
# carried by mult.
# SENSITIVITY: shedding held constant across the symptomatic period.
# This is the limiting case in which viral kinetics contribute nothing and the
# wastewater signal is an infectiousness-weighted count of infectious people.
# Compare against the measured profile in fit_combined.R, which declines
# approximately 3300-fold from stage 1 to stage 3.
shed.I1 <- 10^7.303
shed.I2 <- 10^7.303
shed.I3 <- 10^7.303

shed.I1a <- shed.I1
shed.I2a <- shed.I2
shed.I3a <- shed.I3

shed.I1b <- shed.I1
shed.I2b <- shed.I2
shed.I3b <- shed.I3

shed.A1 <- shed.I1 * tau
shed.A2 <- shed.I2 * tau
shed.A3 <- shed.I3 * tau

shed.A1a <- shed.A1
shed.A2a <- shed.A2
shed.A3a <- shed.A3

shed.A1b <- shed.A1
shed.A2b <- shed.A2
shed.A3b <- shed.A3

alpha_log <- shed.I1 * tau


# ---- Wastewater observations and flow data ----
ww_std <- ww_dat %>% select(log10_cp_per_person_per_day)
ww_obs <- as.numeric(unlist(ww_std$log10_cp_per_person_per_day))
T_WWobs <- length(ww_obs)

flow_L_daily <- ww_dat %>% select(agg_flow_dat)
flow_L_daily <- as.numeric(unlist(flow_L_daily))

# convert flow from litres to millilitres
flow_mL_daily <- flow_L_daily * 1e3

#WW all observation days

flow_all<-ww_flow %>%select(sumflow)
flow_alldaily<-as.numeric(unlist(flow_all))
flow_mlalldaily<-flow_alldaily*1e3


# ---- Wastewater sampling days ----
# observation days are recorded relative to the observation period,
# so shift by burn-in to align with the simulation timeline
raw_ww_days <- ww_dat$Day
ww_sample_days <- raw_ww_days + burn_in_timesteps

# ---- JAGS data list ----
data_list_combined <- list(
  # population and time
  N = N,
  T = T,
  burn_in_timesteps = burn_in_timesteps,
  
  # contact structure
  contact = contact,
  norm_contact_dy = norm_contact_dy,
  
  # initial conditions: unvaccinated
  A2_start = A2_start,
  A3_start = A3_start,
  I2_start = I2_start,
  I3_start = I3_start,
  R_start = R_start,
  Cuminc_start = Cuminc_start,
  
  # initial conditions: first-dose vaccinated
  Va_start = Va_start,
  Pva_start = Pva_start,
  Eva_start = Eva_start,
  A1va_start = A1va_start,
  A2va_start = A2va_start,
  A3va_start = A3va_start,
  I1va_start = I1va_start,
  I2va_start = I2va_start,
  I3va_start = I3va_start,
  
  # initial conditions: second-dose vaccinated
  Vb_start = Vb_start,
  Pvb_start = Pvb_start,
  Evb_start = Evb_start,
  A1vb_start = A1vb_start,
  A2vb_start = A2vb_start,
  A3vb_start = A3vb_start,
  I1vb_start = I1vb_start,
  I2vb_start = I2vb_start,
  I3vb_start = I3vb_start,
  
  # observed data dimensions
  T_WWobs = T_WWobs,
  T_caseobs = T_caseobs,
  
  # vaccination inputs
  alpha_v = alpha_v,
  mu_v = mu_v,
  va_delay = 14,
  vb_delay = 14,
  mu_fixed =0.143,
  
  # shedding inputs
  tau = tau,
  shed_I1 = shed.I1,
  shed_I2 = shed.I2,
  shed_I3 = shed.I3,
  shed_I1a = shed.I1a,
  shed_I2a = shed.I2a,
  shed_I3a = shed.I3a,
  shed_I1b = shed.I1b,
  shed_I2b = shed.I2b,
  shed_I3b = shed.I3b,
  shed_A1 = shed.A1,
  shed_A2 = shed.A2,
  shed_A3 = shed.A3,
  shed_A1a = shed.A1a,
  shed_A2a = shed.A2a,
  shed_A3a = shed.A3a,
  shed_A1b = shed.A1b,
  shed_A2b = shed.A2b,
  shed_A3b = shed.A3b,
  shed_alpha = alpha_log,
  
  # wastewater model inputs
  flow_mL_daily = flow_mL_daily,
  flow_mlalldaily=flow_mlalldaily,
  wwtp_population = 1278020,
  ww_sample_days = ww_sample_days,
  tmax = 15,
  
  # observed data
  ww_obs = ww_obs,
  cases_obs = cases_obsb
)

# ---- Initial values for MCMC chains ----
# Dispersed initial values for two chains
inits_list <- list(
  list(
    beta = 0.8,
    kappa = 0.95,
    report_frac = 0.50,
    log_mult = log(1e-4),
    tau_ww = 0.4,
    transit_time_mean = 2.5,
    transit_time_cv = 0.3,
    .RNG.name = "base::Wichmann-Hill",
    .RNG.seed = 42
  ),
  list(
    beta = 0.9,
    kappa = 0.92,
    report_frac = 0.55,
    log_mult = log(1e-3),
    tau_ww = 0.5,
    transit_time_mean = 2.8,
    transit_time_cv = 0.35,
    .RNG.name = "base::Wichmann-Hill",
    .RNG.seed = 99
  )
)


# ---- Fit combined model ----
system.time({
  combined_fit <- run.jags(
    textstring,
    data = data_list_combined,
    monitor = c(
      "ww_pred", "cases_pred", "log10_conc_all",
      "cases_pred_g", "mu_nb_g",
      "log10_conc_all_g_pred", "log10_conc_all_g",
      "E0", "P0", "A10", "I10",
      "log10_conc", "mu_nb",
      "P_total", "A_total", "I_total",
      "shed_P", "shed_A", "shed_I", "mult", "log_mult",
      "tau_ww", "transit_time_mean", "transit_time_cv",
      "beta", "kappa", "phi",
      "total_Cuminc", "Cuminc", "active_infected", "total_lambda",
      "report_frac", "Vea", "Veb", "m", "sympt_new", "R_total",
      "delta_inv", "theta_invall", "omega_invall", "delayed_conc",
      "daily_shedding", "infectious_pop", "shed_per_inf", "total_new_cases"
    ),
    method = "parallel",
    sample = 40000,
    adapt = 4000,
    burnin = 4000,
    thin = 2,
    n.chains = 2,
    inits = inits_list,
    summarise = FALSE
  )
})

# convert to mcmc.list for downstream diagnostics and plotting
combined_fit_base <- as.mcmc.list(combined_fit)

# save fitted chains
save(combined_fit_base, file = file.path("outputs", "fits", "combined_fit_flatshed_refit.RData"))














