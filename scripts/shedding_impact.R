# ==================================================
# Script: 05_make_manuscript_figures.R
# Purpose: Generate manuscript figures from fitted models
# Project: Mpox_clean
# Updated: 2026-04-22
# ==================================================

rm(list = ls())
set.seed(123)

# ---- Packages ----
library(tidyverse)
library(readxl)
library(coda)
library(bayesplot)
library(ggridges)
library(patchwork)
library(scales)
library(ggpubr)

options(scipen = 999)

# ---- Load fitted model outputs ----
load(file.path("outputs", "fits", "combined_fit_base.RData"))
load(file.path("outputs", "fits", "combined_fit_shed50.RData"))
load(file.path("outputs", "fits", "combined_fit_shed20.RData"))

# ---- Load observed case data ----
case_dat <- read_excel(
  file.path("data", "raw", "case_data_V2.xlsx"),
  sheet = "cases"
)

cases_obsb <- case_dat$total_cases

# ==================================================
# Generate model summaries
# ==================================================

# Convert posterior samples to matrices
mcmc_matrixallbase  <- as.matrix(combined_fit_base)
mcmc_matrixalllo <- as.matrix(combined_fit_shed20)
mcmc_matrixallhi <- as.matrix(combined_fit_shed50)


# ==================================================
# Under-reporting factor (UF) estimation
# ==================================================

compute_UF <- function(mcmc_matrix, cum_reported_cases, burn_in = 30) {
  
  post <- as.matrix(mcmc_matrix)
  
  # Symptomatic incidence
  sympt_cols <- grep("^sympt_new", colnames(post))
  sympt_mat  <- post[, sympt_cols]
  sympt_mat_clean <- sympt_mat[, -c(1:burn_in)]
  total_sympt <- rowSums(sympt_mat_clean, na.rm = TRUE)
  
  UF_sympt <- total_sympt / cum_reported_cases
  
  # Total infections via R compartment
  R_cols <- grep("^R_total\\[", colnames(post))
  R_mat  <- post[, R_cols]
  R_mat_clean <- R_mat[, -c(1:burn_in)]
  R_total_clean <- R_mat_clean[, ncol(R_mat_clean)]
  
  UF_total <- R_total_clean / cum_reported_cases
  
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

# Load observed case data
cum_reported_cases <- sum(cases_obsb)

UF_base <- compute_UF(mcmc_matrixallbase, cum_reported_cases)
UF_lo   <- compute_UF(mcmc_matrixalllo,   cum_reported_cases)
UF_hi   <- compute_UF(mcmc_matrixallhi,   cum_reported_cases)

UF_base$UF_sympt
UF_lo$UF_sympt
UF_hi$UF_sympt

UF_base$UF_total
UF_lo$UF_total
UF_hi$UF_total




