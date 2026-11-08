# ==================================================
# Script: decay_rateimpact.R
# Purpose: Figures and summaries for the viral decay rate sensitivity analysis
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
load(file.path("outputs", "fits", "combined_decay_high.RData"))
load(file.path("outputs", "fits", "combined_decay_low.RData"))

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
mcmc_matrixalllo <- as.matrix(combined_decay_high)
mcmc_matrixallhi <- as.matrix(combined_decay_low)


# E0 = 8 (3–14); P0 = 2 (0–5); A10 = 1 (0–3); I10 = 2 (0–6)


# ---- Posterior summary function ----
# Computes median and 95% credible interval
summary_median_CI <- function(samples) {
  med   <- apply(samples, 2, median)
  lower <- apply(samples, 2, quantile, probs = 0.025)
  upper <- apply(samples, 2, quantile, probs = 0.975)
  
  cbind(
    median       = med,
    lower_95_CI  = lower,
    upper_95_CI  = upper
  )
}


posterior_summarybase <- summary_median_CI(mcmc_matrixallbase)
posterior_summarylo <- summary_median_CI(mcmc_matrixalllo)
posterior_summaryhi <- summary_median_CI(mcmc_matrixallhi)


prev_basesummary <- as.data.frame(posterior_summarybase[grep("active_infected", rownames(posterior_summarybase)), ])
prev_losummary <- as.data.frame(posterior_summarylo[grep("active_infected", rownames(posterior_summarylo)), ])
prev_hisummary <- as.data.frame(posterior_summaryhi[grep("active_infected", rownames(posterior_summaryhi)), ])

burn_in_timesteps <- 30
prev_basesummaryb  <- prev_basesummary[(burn_in_timesteps + 1):nrow(prev_basesummary), ]
prev_losummaryb <- prev_losummary[(burn_in_timesteps + 1):nrow(prev_losummary), ]
prev_hisummaryb <- prev_hisummary[(burn_in_timesteps + 1):nrow(prev_hisummary), ]


plot_baseprev <- data.frame(
  time = 1:nrow(case_dat),
  Date=case_dat$Date,# time index
  median_fit = prev_basesummaryb$median,          # model median fit
  lower_ci = prev_losummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = prev_hisummaryb$upper_95_CI         # upper 95% CI
)


plot_loprev <- data.frame(
  time = 1:nrow(case_dat), 
  Date=case_dat$Date,# time index
  median_fit = prev_losummaryb$median,          # model median fit
  lower_ci = prev_losummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = prev_losummaryb$upper_95_CI         # upper 95% CI
)


plot_hiprev <- data.frame(
  time = 1:nrow(case_dat), 
  Date=case_dat$Date,# time index
  median_fit = prev_hisummaryb$median,          # model median fit
  lower_ci = prev_hisummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = prev_hisummaryb$upper_95_CI         # upper 95% CI
)

#Generate prevalence plots
# Add model labels
plot_baseprev$model <- "base_decay"
plot_loprev$model <- "lo_decay"
plot_hiprev$model <- "hi_decay"

# Combine into one dataframe
plot_allprev <- bind_rows(plot_baseprev, plot_loprev, plot_hiprev)

custom_theme <- theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, colour = "black"),
    axis.title = element_text(size = 14, face = "bold", colour = "black"),
    axis.text = element_text(size = 11, colour = "black")
  )

# plot_allCuminc$model
plot_prev_all <- ggplot(plot_allprev, aes(x = Date, y = median_fit, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +
  geom_line(size = 1.1) +
  labs(
    x = "Date",
    y = "Prevalence(Estimated infected population)",
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("base_decay" = "blue", "lo_decay" = "darkgreen", "hi_decay" = "Orange")) +
  scale_fill_manual(values = c("base_decay" = "blue", "lo_decay" = "darkgreen", "hi_decay" = "Orange")) +
  custom_theme

#Generate incidence plots#
Cuminc_basesummary <- as.data.frame(posterior_summarybase[grep("total_Cuminc", rownames(posterior_summarybase)), ])
Cuminc_losummary <- as.data.frame(posterior_summarylo[grep("total_Cuminc", rownames(posterior_summarylo)), ])
Cuminc_hisummary <- as.data.frame(posterior_summaryhi[grep("total_Cuminc", rownames(posterior_summaryhi)), ])

burn_in_timesteps <- 30
Cuminc_basesummaryb <- Cuminc_basesummary[(burn_in_timesteps + 1):nrow(Cuminc_basesummary), ]
Cuminc_losummaryb <- Cuminc_losummary[(burn_in_timesteps + 1):nrow(Cuminc_losummary), ]
Cuminc_hisummaryb <- Cuminc_hisummary[(burn_in_timesteps + 1):nrow(Cuminc_hisummary), ]


plot_baseCuminc <- data.frame(
  time = 1:nrow(case_dat),
  Date=case_dat$Date,# time index# time index
  median_fit = Cuminc_basesummaryb$median,          # model median fit
  lower_ci = Cuminc_basesummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = Cuminc_basesummaryb$upper_95_CI         # upper 95% CI
)


plot_loCuminc <- data.frame(
  time = 1:nrow(case_dat), 
  Date=case_dat$Date,# time index# time index
  median_fit = Cuminc_losummaryb$median,          # model median fit
  lower_ci = Cuminc_losummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = Cuminc_losummaryb$upper_95_CI         # upper 95% CI
)

plot_hiCuminc <- data.frame(
  time = 1:nrow(case_dat),
  Date=case_dat$Date,# time index# time index
  median_fit = Cuminc_hisummaryb$median,          # model median fit
  lower_ci = Cuminc_hisummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = Cuminc_hisummaryb$upper_95_CI         # upper 95% CI
)


#########generate plots for cumulative incidence
plot_baseCuminc$model <- "base_decay"
plot_loCuminc$model <- "lo_decay"
plot_hiCuminc$model <- "hi_decay"

# Combine into one dataframe
plot_allCuminc <- bind_rows(plot_baseCuminc, plot_loCuminc,plot_hiCuminc)

plot_Cuminc_all <- ggplot(plot_allCuminc, aes(x = Date, y = median_fit, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +
  geom_line(size = 1.1) +
  labs(
    x = "Date",
    y = "Cumulative Incidence",
    #title = "Model-predicted cumulative incidence by model type"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("base_decay" = "blue", "lo_decay" = "darkgreen", "hi_decay" = "Orange")) +
  scale_fill_manual(values = c("base_decay" = "blue", "lo_decay" = "darkgreen", "hi_decay" = "Orange")) +
  custom_theme


decayrate_sens <- plot_prev_all+plot_Cuminc_all +
  #plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


ggsave(
  filename = "outputs/figures/supp_decay_sensitivity.pdf",
  plot = decayrate_sens,
  width = 12,       # adjust width for your manuscript column/page size
  height = 7,      # adjust height accordingly
  units = "in",
  device = cairo_pdf # ensures high-quality text rendering
)


# ==================================================
# Compare Combined (base, low decay, high decay)
# ==================================================

# ------------------------------
# 1. Peak Prevalence
# ------------------------------

burn_in <- 30

# Identify active infected columns
active_cols_base <- grep("active_infected", colnames(mcmc_matrixallbase))
active_cols_lo   <- grep("active_infected", colnames(mcmc_matrixalllo))
active_cols_hi   <- grep("active_infected", colnames(mcmc_matrixallhi))

# Remove burn-in timesteps
active_cols_base <- active_cols_base[-(1:burn_in)]
active_cols_lo   <- active_cols_lo[-(1:burn_in)]
active_cols_hi   <- active_cols_hi[-(1:burn_in)]

# Compute peak per posterior draw
peak_base <- apply(mcmc_matrixallbase[, active_cols_base], 1, max)
peak_lo   <- apply(mcmc_matrixalllo[, active_cols_lo], 1, max)
peak_hi   <- apply(mcmc_matrixallhi[, active_cols_hi], 1, max)

# Posterior probabilities
prob_base_lower_lo  <- mean(peak_base < peak_lo)
prob_base_lower_hi  <- mean(peak_base < peak_hi)

cat("Posterior probability base < low decay:", round(prob_base_lower_lo,3), "\n")
cat("Posterior probability base < high decay:", round(prob_base_lower_hi,3), "\n")

# ------------------------------
# 2. Cumulative Incidence
# ------------------------------

cum_base <- mcmc_matrixallbase[, grep("total_Cuminc", colnames(mcmc_matrixallbase))]
cum_lo   <- mcmc_matrixalllo[,   grep("total_Cuminc", colnames(mcmc_matrixalllo))]
cum_hi   <- mcmc_matrixallhi[,   grep("total_Cuminc", colnames(mcmc_matrixallhi))]

# Ensure equal length
n <- min(length(cum_base), length(cum_lo), length(cum_hi))

cum_base <- cum_base[1:n]
cum_lo   <- cum_lo[1:n]
cum_hi   <- cum_hi[1:n]

# Posterior probabilities
prob_base_lower_lo  <- mean(cum_base < cum_lo)
prob_base_lower_hi  <- mean(cum_base < cum_hi)

cat("Posterior probability base < low decay:", round(prob_base_lower_lo,3), "\n")
cat("Posterior probability base < high decay:", round(prob_base_lower_hi,3), "\n")


cat("Peak prevalence medians:\n")
cat("Base:", median(peak_base), "\n")
cat("Low :", median(peak_lo), "\n")
cat("High:", median(peak_hi), "\n")

cat("Cumulative incidence medians:\n")
cat("Base:", median(cum_base), "\n")
cat("Low :", median(cum_lo), "\n")
cat("High:", median(cum_hi), "\n")


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

# ==================================================
# Transit time sensitivity
# ==================================================

extract_transit <- function(mcmc_obj, label) {
  df <- data.frame(
    value = as.matrix(mcmc_obj)[, "transit_time_mean"],
    model = label
  )
  return(df)
}

transit_df <- bind_rows(
  extract_transit(mcmc_matrixallbase, "μ = 0.143"),
  extract_transit(mcmc_matrixalllo,   "μ = 0.10"),
  extract_transit(mcmc_matrixallhi,   "μ = 0.20")
)

ggplot(transit_df, aes(x = value, fill = model)) +
  geom_density(alpha = 0.4) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Mean wastewater transit time (days)",
    y = "Density"
  ) +
  scale_fill_manual(values = c("steelblue", "darkgreen", "orange"))

compute_transit_summary <- function(mcmc_matrix) {
  transit <- as.matrix(mcmc_matrix)[, "transit_time_mean"]
  c(
    median   = median(transit),
    lower_95 = quantile(transit, 0.025),
    upper_95 = quantile(transit, 0.975)
  )
}

transit_base <- compute_transit_summary(mcmc_matrixallbase)
transit_lo   <- compute_transit_summary(mcmc_matrixalllo)
transit_hi   <- compute_transit_summary(mcmc_matrixallhi)

transit_base
transit_lo
transit_hi


# Table S2 baseline object
median(as.matrix(mcmc_matrixallbase)[, "transit_time_mean"])
quantile(
  as.matrix(mcmc_matrixallbase)[, "transit_time_mean"],
  c(0.025, 0.5, 0.975)
)
















