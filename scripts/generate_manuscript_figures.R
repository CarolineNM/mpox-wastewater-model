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
library(truncnorm)

options(scipen = 999)


# ---- Load observed case data ----
case_dat <- read_excel(
  file.path("data", "raw", "case_data_V2.xlsx"),
  sheet = "cases"
)

cases_obsb <- case_dat$total_cases


# ---- Load observed wastewater data ----
ww_dat <- read_excel(
  file.path("data", "raw", "case_data_V2.xlsx"),
  sheet = "dailyWW"
)

ww_std <- ww_dat %>% select(log10_cp_per_person_per_day)
ww_obs <- as.numeric(unlist(ww_std$log10_cp_per_person_per_day))


# ---- Load fitted model outputs ----
# Refit version: every section reads the refitted key draws below. No .RData
# fits are loaded.

# ==================================================
# Generate model summaries
# ==================================================
# Refitted draws, already matrices
mcmc_matrixallWW <- readRDS(file.path("outputs", "summaries", "refitted_ww_key_draws.rds"))
mcmc_matrixallcas <- readRDS(file.path("outputs", "summaries", "refitted_case_key_draws.rds"))
mcmc_matrixallcom <- readRDS(file.path("outputs", "summaries", "refitted_combined_key_draws.rds"))

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

posterior_summarycas <- summary_median_CI(mcmc_matrixallcas)
posterior_summarycom <- summary_median_CI(mcmc_matrixallcom)
posterior_summaryww  <- summary_median_CI(mcmc_matrixallWW)


# ==================================================
# Extract model-specific summaries
# ==================================================
# ---- Case-only model ----
total_cas_summary      <- as.data.frame(posterior_summarycas[grep("cases_pred", rownames(posterior_summarycas)), ])
total_casb_summary     <- as.data.frame(posterior_summarycas[grep("mu_nb", rownames(posterior_summarycas)), ])
total_cas_ww_summary   <- as.data.frame(posterior_summarycas[grep("^log10_conc\\[", rownames(posterior_summarycas)), ])
total_cas_wwb_summary  <- as.data.frame(posterior_summarycas[grep("ww_pred", rownames(posterior_summarycas)), ])


# ---- WW-only model ----
total_ww_summary       <- as.data.frame(posterior_summaryww[grep("ww_pred", rownames(posterior_summaryww)), ])
total_ww_summaryb      <- as.data.frame(posterior_summaryww[grep("^log10_conc\\[", rownames(posterior_summaryww)), ])
total_ww_cas_summary   <- as.data.frame(posterior_summaryww[grep("mu_nb", rownames(posterior_summaryww)), ])
total_ww_casb_summary  <- as.data.frame(posterior_summaryww[grep("cases_pred", rownames(posterior_summaryww)), ])


# ---- Combined model ----
total_com_ww_summary    <- as.data.frame(posterior_summarycom[grep("ww_pred", rownames(posterior_summarycom)), ])
total_com_ww_summaryb   <- as.data.frame(posterior_summarycom[grep("^log10_conc\\[", rownames(posterior_summarycom)), ])
total_com_cas_summary   <- as.data.frame(posterior_summarycom[grep("^cases_pred\\[", rownames(posterior_summarycom)), ])
total_com_cas_summaryb  <- as.data.frame(posterior_summarycom[grep("^mu_nb\\[", rownames(posterior_summarycom)), ])


# ==================================================
# Generate plotting data
# ==================================================

# ---- WW-only model: WW fit ----
plot_dataww <- data.frame(
  time       = 1:nrow(ww_std),
  Date       = ww_dat$date,
  observed   = ww_obs,
  median_fit = total_ww_summary$median,
  lower_ci   = total_ww_summary$lower_95_CI,
  upper_ci   = total_ww_summary$upper_95_CI
)

plot_datawwb <- data.frame(
  time       = 1:nrow(ww_std),
  Date       = ww_dat$date,
  observed   = ww_obs,
  median_fit = total_ww_summaryb$median,
  lower_ci   = total_ww_summaryb$lower_95_CI,
  upper_ci   = total_ww_summaryb$upper_95_CI
)


# ---- WW-only model: Case predictions ----
plot_dataww_cases <- data.frame(
  time       = 1:nrow(case_dat),
  Date       = case_dat$Date,
  observed   = cases_obsb,
  median_fit = total_ww_cas_summary$median,
  lower_ci   = total_ww_cas_summary$lower_95_CI,
  upper_ci   = total_ww_cas_summary$upper_95_CI
)

plot_dataww_casesb <- data.frame(
  time       = 1:nrow(case_dat),
  Date       = case_dat$Date,
  observed   = cases_obsb,
  median_fit = total_ww_casb_summary$median,
  lower_ci   = total_ww_casb_summary$lower_95_CI,
  upper_ci   = total_ww_casb_summary$upper_95_CI
)


# ---- Case-only model ----
plot_datcas <- data.frame(
  time       = 1:nrow(case_dat),
  Date       = case_dat$Date,
  observed   = cases_obsb,
  median_fit = total_cas_summary$median,
  lower_ci   = total_cas_summary$lower_95_CI,
  upper_ci   = total_cas_summary$upper_95_CI
)

plot_datcasb <- data.frame(
  time       = 1:nrow(case_dat),
  Date       = case_dat$Date,
  observed   = cases_obsb,
  median_fit = total_casb_summary$median,
  lower_ci   = total_casb_summary$lower_95_CI,
  upper_ci   = total_casb_summary$upper_95_CI
)

plot_datacas_ww <- data.frame(
  time       = 1:nrow(ww_std),
  Date       = ww_dat$date,
  observed   = ww_obs,
  median_fit = total_cas_ww_summary$median,
  lower_ci   = total_cas_ww_summary$lower_95_CI,
  upper_ci   = total_cas_ww_summary$upper_95_CI
)

plot_datacas_wwb <- data.frame(
  time       = 1:nrow(ww_std),
  Date       = ww_dat$date,
  observed   = ww_obs,
  median_fit = total_cas_wwb_summary$median,
  lower_ci   = total_cas_wwb_summary$lower_95_CI,
  upper_ci   = total_cas_wwb_summary$upper_95_CI
)


# ---- Combined model ----
plot_datcomcas <- data.frame(
  time       = 1:nrow(case_dat),
  Date       = case_dat$Date,
  observed   = cases_obsb,
  median_fit = total_com_cas_summary$median,
  lower_ci   = total_com_cas_summary$lower_95_CI,
  upper_ci   = total_com_cas_summary$upper_95_CI
)

plot_datcomcasb <- data.frame(
  time       = 1:nrow(case_dat),
  Date       = case_dat$Date,
  observed   = cases_obsb,
  median_fit = total_com_cas_summaryb$median,
  lower_ci   = total_com_cas_summaryb$lower_95_CI,
  upper_ci   = total_com_cas_summaryb$upper_95_CI
)

plot_datcomww <- data.frame(
  time       = 1:nrow(ww_std),
  Date       = ww_dat$date,
  observed   = ww_obs,
  median_fit = total_com_ww_summary$median,
  lower_ci   = total_com_ww_summary$lower_95_CI,
  upper_ci   = total_com_ww_summary$upper_95_CI
)

plot_datcomwwb <- data.frame(
  time       = 1:nrow(ww_std),
  Date       = ww_dat$date,
  observed   = ww_obs,
  median_fit = total_com_ww_summaryb$median,
  lower_ci   = total_com_ww_summaryb$lower_95_CI,
  upper_ci   = total_com_ww_summaryb$upper_95_CI
)

#################Generate posterior plots##########
# ==================================================
# Plot theme and shared geoms
# ==================================================

custom_theme <- theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, colour = "black"),
    axis.title = element_text(size = 14, face = "bold", colour = "black"),
    axis.text  = element_text(size = 11, colour = "black")
  )

# Shared plot layers
plot_geom <- list(
  geom_point(aes(y = observed), colour = "black", size = 1.5, alpha = 0.8),
  geom_line(aes(y = median_fit), colour = "blue4", size = 1.2),
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              fill = "skyblue", alpha = 0.4)
)


# ==================================================
# Case fits
# ==================================================

# Case-only model
plot_casefit <- ggplot(plot_datcas, aes(x = Date)) +
  plot_geom +
  ylim(0, 20) +
  labs(x = "Date", y = "Fit vs. observed cases",
       title = "(Case-only model)") +
  custom_theme

plot_casefitb <- ggplot(plot_datcasb, aes(x = Date)) +
  plot_geom +
  ylim(0, 11) +
  labs(x = "Date", y = "Fit vs. observed cases",
       title = "(Case-only model)") +
  custom_theme


# WW-only model (cases)
plot_casewwfit <- ggplot(plot_dataww_cases, aes(x = Date)) +
  plot_geom +
  ylim(0, 11) +
  labs(x = "Date", y = "Fit vs. observed cases",
       title = "(WW-only model)") +
  custom_theme

plot_casewwfitb <- ggplot(plot_dataww_casesb, aes(x = Date)) +
  plot_geom +
  ylim(0, 20) +
  labs(x = "Date", y = "Fit vs. observed cases",
       title = "(WW-only model)") +
  custom_theme


# Combined model (cases)
plot_comcasefit <- ggplot(plot_datcomcas, aes(x = Date)) +
  plot_geom +
  ylim(0, 20) +
  labs(x = "Date", y = "Fit vs. observed cases",
       title = "(Combined model)") +
  custom_theme

plot_comcasefitb <- ggplot(plot_datcomcasb, aes(x = Date)) +
  plot_geom +
  ylim(0, 11) +
  labs(x = "Date", y = "Fit vs. observed cases",
       title = "(Combined model)") +
  custom_theme


# ==================================================
# Viral load fits
# ==================================================

# WW-only model
plot_wwfit <- ggplot(plot_dataww, aes(x = Date)) +
  plot_geom +
  ylim(0, 8) +
  labs(x = "Date", y = "Fit vs. observed viral load",
       title = "(WW-only model)") +
  custom_theme

plot_wwfitb <- ggplot(plot_datawwb, aes(x = Date)) +
  plot_geom +
  ylim(0, 8) +
  labs(x = "Date", y = "Fit vs. observed viral load",
       title = "(WW-only model)") +
  custom_theme


# Case-only model (viral load)
plot_wwcasefit <- ggplot(plot_datacas_ww, aes(x = Date)) +
  plot_geom +
  coord_cartesian(ylim = c(NA, 8)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8)) +
  labs(x = "Date", y = "Fit vs. observed viral load",
       title = "(Case-only model)") +
  custom_theme

plot_wwcasefitb <- ggplot(plot_datacas_wwb, aes(x = Date)) +
  plot_geom +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8)) +
  labs(x = "Date", y = "Fit vs. observed viral load",
       title = "(Case-only model)") +
  custom_theme


# Combined model (viral load)
plot_comwwfit <- ggplot(plot_datcomww, aes(x = Date)) +
  plot_geom +
  ylim(0, 8) +
  labs(x = "Date", y = "Fit vs. observed viral load",
       title = "(Combined model)") +
  custom_theme

plot_comwwfitb <- ggplot(plot_datcomwwb, aes(x = Date)) +
  plot_geom +
  ylim(0, 8) +
  labs(x = "Date", y = "Fit vs. observed viral load",
       title = "(Combined model)") +
  custom_theme


# ==================================================
# Layout combinations
# ==================================================

ZZ <- (plot_casefit + plot_casewwfitb + plot_comcasefit) /
  (plot_casefitb + plot_casewwfit + plot_comcasefitb)

mm <- ZZ +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


yy <- (plot_wwfit + plot_wwcasefitb + plot_comwwfit) /
  (plot_wwfitb + plot_wwcasefit + plot_comwwfitb)


# Final manuscript figure layout
top_row <- plot_casefitb + plot_casewwfit + plot_comcasefitb
bottom_row <- plot_wwcasefit + plot_wwfitb + plot_comwwfitb

mainfigure1 <- (top_row / bottom_row) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


# Make top row: Case-only, Combined, WW-only
top_rowb <- plot_casefit  + plot_casewwfitb+ plot_comcasefit
bottom_rowb <- plot_wwcasefitb  + plot_wwfit+ plot_comwwfit

# Combine with patchwork, keeping consistent layout
figures1 <- (top_rowb / bottom_rowb) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = "outputs/figures/refit_revfigures1.pdf",
  plot = mainfigure1,
  width = 12,       
  height = 7,      
  units = "in",
  device = cairo_pdf 
)
#supp figure#
ggsave(
  filename = "outputs/figures/refit_revsuppfigures1.pdf",
  plot = figures1,
  width = 12,       
  height = 7,      
  units = "in",
  device = cairo_pdf 
)

ggsave(
  filename = "outputs/figures/refit_suppfigure1.tiff",
  plot = figures1,  # use your actual plot variable here
  width = 12,
  height = 7,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

ggsave(
  filename = "outputs/figures/refit_revfigures1.tiff",
  plot = mainfigure1,  # use your actual plot variable here
  width = 12,
  height = 7,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

#main figure 2#
#for these we use the fig2 fits with tigher Crs#
mcmc_matrixallWWfig2  <- readRDS(file.path("outputs", "summaries", "refitted_ww_key_draws.rds"))
mcmc_matrixallcasfig2 <- readRDS(file.path("outputs", "summaries", "refitted_case_key_draws.rds"))
mcmc_matrixallcomfig2 <- readRDS(file.path("outputs", "summaries", "refitted_combined_key_draws.rds"))

# Function to compute the median and 95% credible interval
summary_median_CI <- function(samples) {
  med <- apply(samples, 2, median)
  lower <- apply(samples, 2, quantile, probs = 0.025)
  upper <- apply(samples, 2, quantile, probs = 0.975)
  summary_table <- cbind(median = med, lower_95_CI = lower, upper_95_CI = upper)
  return(summary_table)
}

posterior_summarycas <- summary_median_CI(mcmc_matrixallcasfig2)
posterior_summarycom <- summary_median_CI(mcmc_matrixallcomfig2)
posterior_summaryww <- summary_median_CI(mcmc_matrixallWWfig2)

prev_wwsummary <- as.data.frame(posterior_summaryww[grep("active_infected", rownames(posterior_summaryww)), ])
prev_cassummary <- as.data.frame(posterior_summarycas[grep("active_infected", rownames(posterior_summarycas)), ])
prev_comsummary <- as.data.frame(posterior_summarycom[grep("active_infected", rownames(posterior_summarycom)), ])

burn_in_timesteps <- 30
prev_wwsummaryb <- prev_wwsummary[(burn_in_timesteps + 1):nrow(prev_wwsummary), ]
prev_cassummaryb <- prev_cassummary[(burn_in_timesteps + 1):nrow(prev_cassummary), ]
prev_comsummaryb <- prev_comsummary[(burn_in_timesteps + 1):nrow(prev_comsummary), ]


plot_wwprev <- data.frame(
  time = 1:nrow(case_dat),
  Date=case_dat$Date,# time index
  median_fit = prev_wwsummaryb$median,          # model median fit
  lower_ci = prev_wwsummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = prev_wwsummaryb$upper_95_CI         # upper 95% CI
)


plot_casprev <- data.frame(
  time = 1:nrow(case_dat), 
  Date=case_dat$Date,# time index
  median_fit = prev_cassummaryb$median,          # model median fit
  lower_ci = prev_cassummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = prev_cassummaryb$upper_95_CI         # upper 95% CI
)


plot_comprev <- data.frame(
  time = 1:nrow(case_dat), 
  Date=case_dat$Date,# time index
  median_fit = prev_comsummaryb$median,          # model median fit
  lower_ci = prev_comsummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = prev_comsummaryb$upper_95_CI         # upper 95% CI
)

#Generate prevalence plots
# Add model labels
plot_wwprev$model <- "WW-only"
plot_casprev$model <- "Case-only"
plot_comprev$model <- "Combined"

# Combine into one dataframe
plot_allprev <- bind_rows(plot_wwprev, plot_casprev, plot_comprev)

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
    y = "Prevalence (Estimated Infected Population)",
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("WW-only" = "blue", "Case-only" = "darkgreen", "Combined" = "Orange")) +
  scale_fill_manual(values = c("WW-only" = "blue", "Case-only" = "darkgreen", "Combined" = "Orange")) +
  custom_theme

#Generate incidence plots#
Cuminc_wwsummary <- as.data.frame(posterior_summaryww[grep("total_Cuminc", rownames(posterior_summaryww)), ])
Cuminc_cassummary <- as.data.frame(posterior_summarycas[grep("total_Cuminc", rownames(posterior_summarycas)), ])
Cuminc_comsummary <- as.data.frame(posterior_summarycom[grep("total_Cuminc", rownames(posterior_summarycom)), ])

burn_in_timesteps <- 30
Cuminc_wwsummaryb <- Cuminc_wwsummary[(burn_in_timesteps + 1):nrow(Cuminc_wwsummary), ]
Cuminc_cassummaryb <- Cuminc_cassummary[(burn_in_timesteps + 1):nrow(Cuminc_cassummary), ]
Cuminc_comsummaryb <- Cuminc_comsummary[(burn_in_timesteps + 1):nrow(Cuminc_comsummary), ]


plot_wwCuminc <- data.frame(
  time = 1:nrow(case_dat),
  Date=case_dat$Date,# time index# time index
  median_fit = Cuminc_wwsummaryb$median,          # model median fit
  lower_ci = Cuminc_wwsummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = Cuminc_wwsummaryb$upper_95_CI         # upper 95% CI
)


plot_casCuminc <- data.frame(
  time = 1:nrow(case_dat), 
  Date=case_dat$Date,# time index# time index
  median_fit = Cuminc_cassummaryb$median,          # model median fit
  lower_ci = Cuminc_cassummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = Cuminc_cassummaryb$upper_95_CI         # upper 95% CI
)

plot_comCuminc <- data.frame(
  time = 1:nrow(case_dat),
  Date=case_dat$Date,# time index# time index
  median_fit = Cuminc_comsummaryb$median,          # model median fit
  lower_ci = Cuminc_comsummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = Cuminc_comsummaryb$upper_95_CI         # upper 95% CI
)


#########generate plots for cumulative incidence
plot_wwCuminc$model <- "WW-only"
plot_casCuminc$model <- "Case-only"
plot_comCuminc$model <- "Combined"

# Combine into one dataframe
plot_allCuminc <- bind_rows(plot_wwCuminc, plot_casCuminc,plot_comCuminc)

plot_Cuminc_all <- ggplot(plot_allCuminc, aes(x = Date, y = median_fit, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +
  geom_line(size = 1.1) +
  labs(
    x = "Date",
    y = "Cumulative Incidence",
    #title = "Model-predicted cumulative incidence by model type"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("WW-only" = "blue", "Case-only" = "darkgreen", "Combined" = "Orange")) +
  scale_fill_manual(values = c("WW-only" = "blue", "Case-only" = "darkgreen", "Combined" = "Orange")) +
  custom_theme


mainfigure2 <- plot_prev_all+plot_Cuminc_all +
  #plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = "outputs/figures/refit_mainfigure2b.pdf",
  plot = mainfigure2,
  width = 12,       # adjust width for your manuscript column/page size
  height = 7,      # adjust height accordingly
  units = "in",
  device = cairo_pdf # ensures high-quality text rendering
)

ggsave(
  filename = "outputs/figures/refit_mainfigure2b.pdf",
  plot = mainfigure2,
  width = 12,       # adjust width for your manuscript column/page size
  height = 7,      # adjust height accordingly
  units = "in",
  device = cairo_pdf # ensures high-quality text rendering
)

ggsave(
  filename = "outputs/figures/refit_mainfigure2b.tiff",
  plot = mainfigure2,  # use your actual plot variable here
  width = 12,
  height = 7,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

# ==================================================
# Posterior probability comparisons
# ==================================================
# Compare prevalence and cumulative incidence across:
#   - WW-only model
#   - Case-only model
#   - Combined model


# ==================================================
# 1. Prevalence (active infected)
# ==================================================

# Extract posterior draws
prev_WW  <- mcmc_matrixallWWfig2[,  grep("active_infected", colnames(mcmc_matrixallWWfig2))]
prev_cas <- mcmc_matrixallcasfig2[, grep("active_infected", colnames(mcmc_matrixallcasfig2))]
prev_com <- mcmc_matrixallcomfig2[, grep("active_infected", colnames(mcmc_matrixallcomfig2))]

# Ensure equal length
n <- min(length(prev_WW), length(prev_cas), length(prev_com))

prev_WW  <- prev_WW[1:n]
prev_cas <- prev_cas[1:n]
prev_com <- prev_com[1:n]

# Posterior probabilities (overall prevalence)
prob_com_lower_WW  <- mean(prev_com < prev_WW)
prob_com_lower_cas <- mean(prev_com < prev_cas)

cat("Posterior probability that combined < wastewater-only:",
    round(prob_com_lower_WW, 3), "\n")

cat("Posterior probability that combined < case-only:",
    round(prob_com_lower_cas, 3), "\n")


# ==================================================
# 2. Peak prevalence (excluding burn-in)
# ==================================================

# Identify active infected columns
active_cols <- grep("active_infected", colnames(mcmc_matrixallcom))

# Remove first 30 timesteps (burn-in)
active_cols_burn <- active_cols[-(1:30)]

# Compute peak prevalence per posterior draw
peakprev_com <- apply(mcmc_matrixallcom[, active_cols_burn], 1, max)
peakprev_WW  <- apply(mcmc_matrixallWW[,  active_cols_burn], 1, max)
peakprev_cas <- apply(mcmc_matrixallcas[, active_cols_burn], 1, max)

# Posterior probabilities (peak prevalence)
prob_com_lower_WW  <- mean(peakprev_com < peakprev_WW)
prob_com_lower_cas <- mean(peakprev_com < peakprev_cas)

cat("Posterior probability that combined < wastewater-only:",
    round(prob_com_lower_WW, 3), "\n")

cat("Posterior probability that combined < case-only:",
    round(prob_com_lower_cas, 3), "\n")


# ==================================================
# 3. Cumulative incidence
# ==================================================

# Extract posterior draws
Cuminc_WW  <- mcmc_matrixallWW[,  grep("total_Cuminc", colnames(mcmc_matrixallWW))]
Cuminc_cas <- mcmc_matrixallcas[, grep("total_Cuminc", colnames(mcmc_matrixallcas))]
Cuminc_com <- mcmc_matrixallcom[, grep("total_Cuminc", colnames(mcmc_matrixallcom))]

# Ensure equal length
n <- min(length(Cuminc_WW), length(Cuminc_cas), length(Cuminc_com))

Cuminc_WW  <- Cuminc_WW[1:n]
Cuminc_cas <- Cuminc_cas[1:n]
Cuminc_com <- Cuminc_com[1:n]

# Posterior probabilities (cumulative incidence)
cumin_com_lower_WW  <- mean(Cuminc_com < Cuminc_WW)
cumin_com_lower_cas <- mean(Cuminc_com < Cuminc_cas)

cat("Posterior probability that combined < wastewater-only:",
    round(cumin_com_lower_WW, 3), "\n")

cat("Posterior probability that combined < case-only:",
    round(cumin_com_lower_cas, 3), "\n")


#main figure 3b#
# Refit version: mcmc_matrixallcom already holds the refitted combined draws
# from the top of the script.

# Panel B needs log10_conc_all_g and mu_nb_g, which are not in the key draws.
# They come from a separate file written by extract_extra_series.R. If it is
# not there yet, panel B is skipped and the rest of the script carries on.
extra_file_com  <- file.path("outputs", "summaries", "refitted_combined_extra_draws.rds")
extra_available <- file.exists(extra_file_com)

if (extra_available) {
  extra_com <- readRDS(extra_file_com)
  if (nrow(extra_com) != nrow(mcmc_matrixallcom)) {
    stop("refitted_combined_extra_draws.rds has ", nrow(extra_com),
         " draws but the key draws have ", nrow(mcmc_matrixallcom),
         ". They must be thinned identically to be combined.")
  }
  mcmc_matrixallcom <- cbind(mcmc_matrixallcom, extra_com)
} else {
  message("NOTE: ", extra_file_com, " not found. ",
          "Figure 3 panel B (sliding correlation) is skipped; ",
          "figure 3 is saved without it.")
}

summary_median_CI <- function(samples) {
  med <- apply(samples, 2, median)
  lower <- apply(samples, 2, quantile, probs = 0.025)
  upper <- apply(samples, 2, quantile, probs = 0.975)
  summary_table <- cbind(median = med, lower_95_CI = lower, upper_95_CI = upper)
  return(summary_table)
}

posterior_summarycom <- summary_median_CI(mcmc_matrixallcom)

# ---- 1) compute correlations for each group ----
plot_sliding_corr_by_group <- function(mcmc_matrix, posterior_summary, groups = 1:3, 
                                       burn_in_timesteps = 30,
                                       ww_prefix = "^log10_conc_all_g\\[",
                                       case_prefix = "^mu_nb_g\\[") {
  
  all_summary <- list()
  
  for (g in groups) {
    # 1. Extract matrices for WW and cases
    ww_mat <- mcmc_matrix[, grep(paste0(ww_prefix, g, ","), colnames(mcmc_matrix))]
    case_mat <- mcmc_matrix[, grep(paste0(case_prefix, g, ","), colnames(mcmc_matrix))]
    
    # 2. Define window sizes
    window_sizes <- seq(7, ncol(ww_mat), by = 7)
    n_draws <- nrow(ww_mat)
    
    # 3. Loop over window sizes
    correlation_list <- vector("list", length(window_sizes))
    names(correlation_list) <- window_sizes
    
    for (w in seq_along(window_sizes)) {
      size <- window_sizes[w]
      cor_per_draw <- numeric(n_draws)
      
      for (d in 1:n_draws) {
        ww_sub <- ww_mat[d, 1:size]
        case_sub <- case_mat[d, 1:size]
        
        if (sd(ww_sub) > 0 && sd(case_sub) > 0) {
          cor_per_draw[d] <- cor(ww_sub, case_sub)
        } else {
          cor_per_draw[d] <- NA
        }
      }
      correlation_list[[w]] <- cor_per_draw
    }
    
    # 4. Label for group
    group_label <- c(
      "Group 1 (low activity)",
      "Group 2 (medium activity)",
      "Group 3 (high activity)"
    )[g]
    
    # 5. Summarise results
    summary_df <- do.call(rbind, lapply(seq_along(correlation_list), function(i) {
      vals <- na.omit(correlation_list[[i]])
      data.frame(
        group =  group_label,
        window_size = window_sizes[i],
        frac_time = window_sizes[i] / max(window_sizes),
        median = median(vals),
        lower = quantile(vals, 0.025),
        upper = quantile(vals, 0.975)
      )
    }))
    
    # 6. Peak transmission fraction
    FOI <- as.data.frame(posterior_summary[grep("^total_lambda\\[", rownames(posterior_summary)), "median"])
    FOI_final <- FOI[(burn_in_timesteps + 1):nrow(FOI), ]
    transmission_peak_day <- which.max(FOI_final)
    total_days <- ncol(ww_mat)
    transmission_peak_frac <- transmission_peak_day / total_days
    
    summary_df$transmission_peak_frac <- transmission_peak_frac
    all_summary[[g]] <- summary_df
  }
  
  # Combine and return
  all_summary_df <- bind_rows(all_summary)
  return(all_summary_df)
}

# Usage:
if (extra_available) {
  all_summary_df <- plot_sliding_corr_by_group(mcmc_matrixallcom, posterior_summarycom)
}

shedP_1 <- mcmc_matrixallcom[, grep("shed_P\\[1,", colnames(mcmc_matrixallcom))]
shedP_2 <- mcmc_matrixallcom[, grep("shed_P\\[2,", colnames(mcmc_matrixallcom))]
shedP_3 <- mcmc_matrixallcom[, grep("shed_P\\[3,", colnames(mcmc_matrixallcom))]

shedA_1 <- mcmc_matrixallcom[, grep("shed_A\\[1,", colnames(mcmc_matrixallcom))]
shedA_2 <- mcmc_matrixallcom[, grep("shed_A\\[2,", colnames(mcmc_matrixallcom))]
shedA_3 <- mcmc_matrixallcom[, grep("shed_A\\[3,", colnames(mcmc_matrixallcom))]

shedI_1 <- mcmc_matrixallcom[, grep("shed_I\\[1,", colnames(mcmc_matrixallcom))]
shedI_2 <- mcmc_matrixallcom[, grep("shed_I\\[2,", colnames(mcmc_matrixallcom))]
shedI_3 <- mcmc_matrixallcom[, grep("shed_I\\[3,", colnames(mcmc_matrixallcom))]


# Total shedding for each group (denominator)

shed1_total <- shedP_1 + shedA_1 + shedI_1
shed2_total <- shedP_2 + shedA_2 + shedI_2
shed3_total <- shedP_3 + shedA_3 + shedI_3
shedall_total <- shed1_total + shed2_total + shed3_total


# Relative contribution of each group
relg_1 <- shed1_total / shedall_total
relg_2 <- shed2_total/ shedall_total
relg_3 <- shed3_total / shedall_total

# Function to summarize posterior samples over time
summary_rel_group <- function(mat, group_label) {
  data.frame(
    time = 1:ncol(mat),
    Date=case_dat$Date,# time index
    median = apply(mat, 2, median, na.rm = TRUE),
    lower = apply(mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(mat, 2, quantile, probs = 0.975, na.rm = TRUE),
    group = group_label
  )
}

burn_in <- 30

# Remove first 30 timesteps from each matrix
relg_1_burned <- relg_1[, (burn_in):ncol(relg_1)]
relg_2_burned <- relg_2[, (burn_in):ncol(relg_2)]
relg_3_burned <- relg_3[, (burn_in):ncol(relg_3)]

# Combine summaries for plotting
relg_all <- dplyr::bind_rows(
  summary_rel_group(relg_1_burned, "Group 1 (low activity)"),
  summary_rel_group(relg_2_burned, "Group 2 (medium activity)"),
  summary_rel_group(relg_3_burned, "Group 3 (high activity)")
)


plotshed=ggplot(relg_all, aes(x = Date, y = median, fill = group, color = group)) +
  geom_line(size = 1) +
  geom_point(size = 0, alpha = 0)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, linetype = 0) +
  labs(
    #title = "Relative contribution to shedding by sexual activity group",
    x = "Time",
    y = "Proportion of viral load shed"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(legend.position = "bottom")

plotshed

#We now do the same for the total number of infected individuals
P1_total <- mcmc_matrixallcom[, grep("P_total\\[1,", colnames(mcmc_matrixallcom))]
P2_total <- mcmc_matrixallcom[, grep("P_total\\[2,", colnames(mcmc_matrixallcom))]
P3_total <- mcmc_matrixallcom[, grep("P_total\\[3,", colnames(mcmc_matrixallcom))]

A1_total <- mcmc_matrixallcom[, grep("A_total\\[1,", colnames(mcmc_matrixallcom))]
A2_total <- mcmc_matrixallcom[, grep("A_total\\[2,", colnames(mcmc_matrixallcom))]
A3_total <- mcmc_matrixallcom[, grep("A_total\\[3,", colnames(mcmc_matrixallcom))]

I1_total <- mcmc_matrixallcom[, grep("I_total\\[1,", colnames(mcmc_matrixallcom))]
I2_total <- mcmc_matrixallcom[, grep("I_total\\[2,", colnames(mcmc_matrixallcom))]
I3_total <- mcmc_matrixallcom[, grep("I_total\\[3,", colnames(mcmc_matrixallcom))]


# Sum across compartments per group
infected_1 <- P1_total + A1_total + I1_total
infected_2 <- P2_total + A2_total + I2_total
infected_3 <- P3_total + A3_total + I3_total

burn_in <- 30  # Number of initial time steps to exclude

# Remove burn-in columns
infected_1_burned <- infected_1[, (burn_in):ncol(infected_1)]
infected_2_burned <- infected_2[, (burn_in):ncol(infected_2)]
infected_3_burned <- infected_3[, (burn_in):ncol(infected_3)]

# Summary function for CrI and median over time
summarise_draws <- function(mat, group_label) {
  data.frame(
    time = 1:ncol(mat),
    Date=case_dat$Date,# time index
    median = apply(mat, 2, median),
    lower = apply(mat, 2, quantile, probs = 0.025),
    upper = apply(mat, 2, quantile, probs = 0.975),
    group = group_label
  )
}
# Apply summarization after trimming
infected_df <- bind_rows(
  summarise_draws(infected_1_burned, "Group 1 (low activity)"),
  summarise_draws(infected_2_burned, "Group 2 (medium activity)"),
  summarise_draws(infected_3_burned, "Group 3 (high activity)")
)

relg_v=infected_df %>% filter(group=="Group 1 (low activity)")
relg_vv=infected_df %>% filter(group=="Group 2 (medium activity)")
relg_vvv=infected_df %>% filter(group=="Group 3 (high activity)")
max(relg_vvv$median)

custom_themeb <- theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, colour = "black"),
    axis.title = element_text(size = 16, face = "bold", colour = "black"),
    legend.text=element_text(size = 16, face = "bold", colour = "black"),
    legend.title=element_text(size = 16, face = "bold", colour = "black"),
    axis.text = element_text(size = 16, colour = "black")
  )

custom_themec <- theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, colour = "black"),
    
    axis.title = element_text(size = 16, face = "bold", colour = "black"),
    axis.text  = element_text(size = 14, colour = "black"),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14),
    
    legend.margin = margin(t = 10, b = 5),
    legend.box.margin = margin(t = 10),
    legend.spacing.x = unit(0.6, "cm"),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey85")
  )

custom_theme <- theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, colour = "black"),
    axis.title = element_text(size = 16, face = "bold", colour = "black"),
    legend.text=element_text(size = 16, face = "bold", colour = "black"),
    legend.title=element_text(size = 16, face = "bold", colour = "black"),
    axis.text = element_text(size = 16, colour = "black")
  )


group_colors <- c(
  "Group 1 (low activity)" = "#e41a1c",
  "Group 2 (medium activity)" = "#4daf4a",
  "Group 3 (high activity)" = "#377eb8"
)
# Plotting
# Plot 1: Relative contribution to viral shedding
plot_shed <- ggplot(relg_all, aes(x = Date, y = median, fill = group, color = group)) +
  geom_line(size = 1.1) +
  geom_point(size = 0, alpha = 0)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  labs(
    x = "Date",
    y = "Proportion of total viral load"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = group_colors, name = NULL) +
  scale_fill_manual(values = group_colors, name = NULL) +
  theme_minimal(base_size = 13) +
  custom_themeb


# Keep legend only in p_corr
# Plot 2: Number of infectious individuals by group
plot_epi <- ggplot(infected_df, aes(x = Date, y = median, color = group, fill = group)) +
  geom_line(size = 1.1) +
  geom_point(size = 0, alpha = 0)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  labs(
    #title = "Number of infectious individuals by sexual activity group",
    x = "Date",
    y = "Estimated number of infectious individuals"
  ) +
  scale_color_manual(values = group_colors, name = NULL) +
  scale_fill_manual(values = group_colors, name = NULL) +
  theme_minimal(base_size = 13) +
  custom_theme

group_colors <- c(
  "Group 1 (low activity)" = "#e41a1c",
  "Group 2 (medium activity)" = "#4daf4a",
  "Group 3 (high activity)" = "#377eb8"
)

plot_sliding_corr <- function(all_summary_df, group_colors) {
  
  p <- ggplot(all_summary_df, aes(x = frac_time, y = median, color = group, fill = group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    geom_line(size = 1) +
    geom_point(size = 1.5, show.legend = FALSE) + # hide point legend
    geom_vline(aes(xintercept = transmission_peak_frac), linetype = "dashed", color = "black") +
    annotate("text", x = 0.27, y = 0.85, label = "Peak transmission", hjust = 0, size = 3) +
    facet_wrap(~ group, scales = "free_y", ncol = 1) +
    scale_color_manual(values = group_colors, name = "group") +
    scale_fill_manual(values = group_colors, name = "group") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = "Fraction of time series used",
      y = "Pearson correlation (viral load vs. cases)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  return(p)
}


plot_sliding_corrb <- function(all_summary_df, group_colors) {
  
  peak_labels <- data.frame(
    group = c("Group 1 (low activity)",
              "Group 2 (medium activity)",
              "Group 3 (high activity)"),
    x = c(0.27, 0.27, 0.27),
    y = c(0.80, 0.70, 0.90),
    label = "Peak transmission"
  )
  
  p <- ggplot(all_summary_df,
              aes(x = frac_time,
                  y = median,
                  color = group,
                  fill = group)) +
    
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.2,
                color = NA) +
    
    geom_line(size = 1) +
    geom_point(size = 1.5, show.legend = FALSE) +
    
    geom_vline(aes(xintercept = transmission_peak_frac),
               linetype = "dashed",
               color = "black") +
    
    geom_text(data = peak_labels,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              size = 3) +
    
    facet_wrap(~ group, scales = "free_y", ncol = 1) +
    
    scale_color_manual(values = group_colors, name = "Group") +
    scale_fill_manual(values = group_colors, name = "Group") +
    
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    
    labs(
      x = "Fraction of time series used",
      y = "Pearson correlation (viral load vs. cases)"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  return(p)
}
if (extra_available) {
  p_corr <- plot_sliding_corrb(all_summary_df, group_colors)
  p_corr  <- p_corr + custom_themeb
}

# Combine with  other two plots:
if (extra_available) {
  top_row <- plot_shed + p_corr + plot_layout(widths = c(1, 1))
} else {
  top_row <- plot_shed
}
mainfigure3 <- top_row / plot_epi +
  plot_annotation(tag_levels = "A")  
mainfigure3

ggsave(
  filename = "outputs/figures/refit_mainfigure3.pdf",
  plot = mainfigure3,
  width = 14,       
  height = 11,      
  units = "in",
  device = cairo_pdf # ensures high-quality text rendering
)

ggsave(
  filename = "outputs/figures/refit_mainfigure3.tiff",
  plot = mainfigure3,  # use your actual plot variable here
  width = 14,
  height = 11,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

#Marginal posterior plots#
case_mod  <- mcmc_matrixallcas
ww_mod   <- mcmc_matrixallWW
com_mod <- mcmc_matrixallcom

###Define parameter to model mapping
param_model_map <- list(
  beta = c("Case-only", "Viral load-only", "Combined"),
  kappa = c("Case-only", "Viral load-only", "Combined"),
  delta_inv=c("Case-only", "Viral load-only", "Combined"),
  theta_invall=c("Case-only", "Viral load-only", "Combined"),
  omega_invall=c("Case-only", "Viral load-only", "Combined"),
  Vea=c("Case-only", "Viral load-only", "Combined"),
  Veb= c("Case-only", "Viral load-only", "Combined"),
  phi = c("Case-only","Viral load-only", "Combined"),
  report_frac = c("Case-only","Viral load-only","Combined"),
  m = c("Case-only","Viral load-only","Combined"),
  transit_time_mean = c("Case-only","Viral load-only", "Combined"),
  transit_time_cv = c("Case-only","Viral load-only", "Combined"),
  transit_time_cv = c("Case-only","Viral load-only", "Combined"),
  mult = c("Case-only","Viral load-only", "Combined"),
  tau_ww = c("Viral load-only","Case-only", "Combined")
)

#Add parameter labels
param_labels <- c(
  beta = "Transmission probability, β",
  kappa = "Mixing probability, κ",
  phi = "Negative binomial dispersion, φ",
  report_frac = "Reporting fraction",
  m = "Proportion asymptomatic, m",
  delta_inv = "Latent period (1/δ), days",
  theta_invall = "Infectious period (1/θ), days",
  omega_invall = "Asymptomatic period (1/ω), days",
  Vea = "Vaccine effectiveness (1st dose)",
  Veb = "Vaccine effectiveness (2nd dose)",
  transit_time_mean = "Mean transit time (days)",
  transit_time_cv = "Sd of transit time",
  mult = "Scaling factor for viral load",
  tau_ww = "Precision of WW likelihood"
)

#Set prior definitions
prior_definitions <- list(
  beta = function(n) rtruncnorm(n, 0, 1, mean = 0.8, sd = 0.1),
  kappa = function(n) rbeta(n, 40, 2),
  phi = function(n) rgamma(n, shape = 2, rate = 0.5),
  report_frac = function(n) rbeta(n, 10, 10),
  m = function(n) rbeta(n, 49.09, 94.88),
  delta_inv=function(n) rgamma(n, shape = 300.6, rate = 47.98),
  theta_invall=function(n) rgamma(n, shape = 42.3, rate = 1.99),
  omega_invall=function(n) rgamma(n, shape = 20.6, rate = 1.45),
  Vea=function(n) rbeta(n, shape1 = 49.3, shape2 = 87.4),
  Veb=function(n) rbeta(n, shape1 = 69.3, shape2 = 35.6),
  transit_time_mean = function(n) rtruncnorm(n, a = 0.1, b = 10,mean = 2.5, sd = sqrt(1 / 0.25)),
  transit_time_cv = function(n) rtruncnorm(n, a = 0.1, b = 1, mean = 0.3, sd = sqrt(1 / 3)),
  mult = function(n) exp(rtruncnorm(n, a = log(1e-8), b = log(1e0), mean = log(1.5e-4), sd = sqrt(1 / 0.05))),
  tau_ww = function(n) rgamma(n, shape = 40, rate = 48)
)


extract_param <- function(mcmc_obj, param_name, model_label) {
  df <- as.data.frame(as.matrix(mcmc_obj)[, param_name])
  colnames(df) <- "value"
  df$source <- model_label
  df$parameter <- param_name
  return(df)
}

#compile full dataset#
n_samples <- 40000
models <- list("Case-only" = case_mod, "Viral load-only" = ww_mod, "Combined" = com_mod)

all_param_dfs <- list()

for (param in names(param_model_map)) {
  # Simulate prior
  prior_df <- data.frame(
    value = prior_definitions[[param]](n_samples),
    source = "Prior",
    parameter = param
  )
  
  # Posteriors
  posteriors <- list()
  for (model_name in param_model_map[[param]]) {
    posteriors[[model_name]] <- extract_param(models[[model_name]], param, model_name)
  }
  
  all_param_dfs[[param]] <- bind_rows(prior_df, do.call(bind_rows, posteriors))
}

beta_df_all <- bind_rows(all_param_dfs)
beta_df_all$source <- factor(beta_df_all$source, levels = c("Combined","Viral load-only","Case-only","Prior"))
beta_df_all$parameter_label <- param_labels[beta_df_all$parameter]

beta_df_allb <- beta_df_all %>%
  mutate(
    value_plot = case_when(
      parameter == "mult" ~ log10(value),
      TRUE ~ value
    ),
    
    parameter_label = case_when(
      parameter == "mult" ~ "log10(ψ)",
      TRUE ~ parameter_label
    ),
    source=case_when(
      source == "Viral load-only" ~ "Wastewater-only",
      TRUE ~ source
    ),
  )
WW_params <- c("mult", "transit_time_mean","tau_ww")

beta_df_allb$source <- factor(beta_df_allb$source, levels = c("Combined","Wastewater-only","Case-only","Prior"))

#plot a,d and c will be based on these fits#
plota=ggplot(beta_df_allb %>% filter(parameter %in% WW_params), 
             aes(x = value_plot, y = source, fill = source)) +
  geom_density_ridges(alpha = 0.7, 
                      quantile_lines = TRUE, 
                      quantiles = 0.5,
                      scale = 0.9,
                      rel_min_height = 0.0005) +
 facet_wrap(
  ~ parameter_label,
  scales = "free_x",
  ncol = 3,
  labeller = labeller(
    parameter_label = c(
      "log[10] (ψ)" = expression(log[10](psi))
    )
  )
)+
  scale_fill_manual(values = c("gray30", "gray60", "gray70", "gray90")) +
  theme_minimal(base_size = 16) +
  ylab("Source")+
  xlab("Parameter value")+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  )
plota

case_params <- c("beta", "kappa","report_frac")

plotbb=ggplot(beta_df_allb %>% filter(parameter %in% case_params), 
             aes(x = value, y = source, fill = source)) +
  geom_density_ridges(alpha = 0.7, 
                      quantile_lines = TRUE, 
                      quantiles = 0.5,
                      scale = 0.9,
                      rel_min_height = 0.0005) +
  facet_wrap(~ parameter_label, scales = "free_x", ncol = 4) +
  scale_fill_manual(values = c("gray30", "gray60", "gray70", "gray90")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  )
plotbb


Inf_params <- c("delta_inv", "theta_invall","omega_invall")
plotc=ggplot(beta_df_all %>% filter(parameter %in% Inf_params), 
             aes(x = value, y = source, fill = source)) +
  geom_density_ridges(alpha = 0.7, 
                      quantile_lines = TRUE, 
                      quantiles = 0.5,
                      scale = 0.9,
                      rel_min_height = 0.0005) +
  facet_wrap(~ parameter_label, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = c("gray30", "gray60", "gray70", "gray90")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  )
plotc


other_params <- c("Vea", "Veb","m")
plotd=ggplot(beta_df_all %>% filter(parameter %in% other_params), 
             aes(x = value, y = source, fill = source)) +
  geom_density_ridges(alpha = 0.7, 
                      quantile_lines = TRUE, 
                      quantiles = 0.5,
                      scale = 0.9,
                      rel_min_height = 0.0005) +
  facet_wrap(~ parameter_label, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = c("gray30", "gray60", "gray70", "gray90")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  )
plotd


#########plot b will be based on fits_fig2#########
case_mod  <- mcmc_matrixallcasfig2
ww_mod   <- mcmc_matrixallWWfig2
com_mod <- mcmc_matrixallcomfig2

###Define parameter to model mapping
param_model_map <- list(
  beta = c("Case-only", "Viral load-only", "Combined"),
  kappa = c("Case-only", "Viral load-only", "Combined"),
  delta_inv=c("Case-only", "Viral load-only", "Combined"),
  theta_invall=c("Case-only", "Viral load-only", "Combined"),
  omega_invall=c("Case-only", "Viral load-only", "Combined"),
  Vea=c("Case-only", "Viral load-only", "Combined"),
  Veb= c("Case-only", "Viral load-only", "Combined"),
  phi = c("Case-only","Viral load-only", "Combined"),
  report_frac = c("Case-only","Viral load-only","Combined"),
  m = c("Case-only","Viral load-only","Combined"),
  transit_time_mean = c("Case-only","Viral load-only", "Combined"),
  transit_time_cv = c("Case-only","Viral load-only", "Combined"),
  transit_time_cv = c("Case-only","Viral load-only", "Combined"),
  mult = c("Case-only","Viral load-only", "Combined"),
  tau_ww = c("Viral load-only","Case-only", "Combined")
)

###########Add parameter labels
param_labels <- c(
  beta = "Transmission probability, β",
  kappa = "Mixing probability, κ",
  phi = "Negative binomial dispersion, φ",
  report_frac = "Reporting fraction",
  m = "Proportion asymptomatic, m",
  delta_inv = "Latent period (1/δ), days",
  theta_invall = "Infectious period (1/θ), days",
  omega_invall = "Asymptomatic period (1/ω), days",
  Vea = "Vaccine effectiveness (1st dose)",
  Veb = "Vaccine effectiveness (2nd dose)",
  transit_time_mean = "Mean transit time (days)",
  transit_time_cv = "Sd of transit time",
  mult = "Scaling factor for viral load",
  tau_ww = "Precision of WW likelihood (τ)"
)


# # # Set prior definitions
prior_definitions <- list(
  beta = function(n) rtruncnorm(n, 0, 1, mean = 0.8, sd = 0.1),
  kappa = function(n) rbeta(n, 40, 2),
  phi = function(n) rgamma(n, shape = 2, rate = 0.5),
  report_frac = function(n) rbeta(n, 10, 10),
  m = function(n) rbeta(n, 49.09, 94.88),
  delta_inv=function(n) rgamma(n, shape = 300.6, rate = 47.98),
  theta_invall=function(n) rgamma(n, shape = 42.3, rate = 1.99),
  omega_invall=function(n) rgamma(n, shape = 20.6, rate = 1.45),
  Vea=function(n) rbeta(n, shape1 = 49.3, shape2 = 87.4),
  Veb=function(n) rbeta(n, shape1 = 69.3, shape2 = 35.6),
  transit_time_mean = function(n) rtruncnorm(n, a = 0.1, b = 10,mean = 2.5, sd = sqrt(1 / 0.25)),
  transit_time_cv = function(n) rtruncnorm(n, a = 0.1, b = 1, mean = 0.3, sd = sqrt(1 / 3)),
  mult = function(n) exp(rtruncnorm(n, a = log(1e-8), b = log(1e0), mean = log(1.5e-4), sd = sqrt(1 / 0.05))),
  tau_ww = function(n) rgamma(n, shape = 40, rate = 48)
)


extract_param <- function(mcmc_obj, param_name, model_label) {
  df <- as.data.frame(as.matrix(mcmc_obj)[, param_name])
  colnames(df) <- "value"
  df$source <- model_label
  df$parameter <- param_name
  return(df)
}

##########compile full dataset
n_samples <- 40000
models <- list("Case-only" = case_mod, "Viral load-only" = ww_mod, "Combined" = com_mod)

all_param_dfs <- list()

for (param in names(param_model_map)) {
  # Simulate prior
  prior_df <- data.frame(
    value = prior_definitions[[param]](n_samples),
    source = "Prior",
    parameter = param
  )
  
  # Posteriors
  posteriors <- list()
  for (model_name in param_model_map[[param]]) {
    posteriors[[model_name]] <- extract_param(models[[model_name]], param, model_name)
  }
  
  all_param_dfs[[param]] <- bind_rows(prior_df, do.call(bind_rows, posteriors))
}

beta_df_all <- bind_rows(all_param_dfs)
beta_df_all$source <- factor(beta_df_all$source, levels = c("Combined","Viral load-only","Case-only","Prior"))
beta_df_all$parameter_label <- param_labels[beta_df_all$parameter]


beta_df_allc <- beta_df_all %>%
  mutate(
    source=case_when(
      source == "Viral load-only" ~ "Wastewater-only",
      TRUE ~ source
    ),
  )

beta_df_allc$source <- factor(beta_df_allc$source, levels = c("Combined","Wastewater-only","Case-only","Prior"))

case_params <- c("beta", "kappa","report_frac")

plotb=ggplot(beta_df_allc %>% filter(parameter %in% case_params), 
             aes(x = value, y = source, fill = source)) +
  geom_density_ridges(alpha = 0.7, 
                      quantile_lines = TRUE, 
                      quantiles = 0.5,
                      scale = 0.9,
                      rel_min_height = 0.0005) +
  facet_wrap(~ parameter_label, scales = "free_x", ncol = 4) +
  scale_fill_manual(values = c("gray30", "gray60", "gray70", "gray90")) +
  theme_minimal(base_size = 16) +
  ylab("Source")+
  xlab("Parameter value")+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  )
plotb
#We then combine all plots#
## Row 1: Reduced Shedding
row1 <- ggpubr::ggarrange(
  plota,
  labels = "B")

row1 <- annotate_figure(row1,
                        top = text_grob("Wastewater parameters", face = "bold", size = 14)
)

row1


row2 <- ggpubr::ggarrange(
  plotb,
  labels = "A")

row2 <- annotate_figure(row2,
                        top = text_grob("Transmission and surveillance parameters", face = "bold", size = 14)
)

row2


row3 <- ggpubr::ggarrange(
  plotc,
  labels = "A")

row3 <- annotate_figure(row3,
                        top = text_grob("Infection parameters", face = "bold", size = 14)
)

row3

row4 <- ggpubr::ggarrange(
  plotd,
  labels = "B")

row4 <- annotate_figure(row4,
                        top = text_grob("Vaccine/natural history parameters", face = "bold", size = 14)
)

row4


finals2 <- ggarrange(
  row2, row1,
  ncol= 1,
  heights = c(1, 1)
)
finals2 

finals3 <- ggarrange(
  row3, row4,
  ncol= 1,
  heights = c(1, 1)
)
finals3 


ggsave(
  filename = "outputs/figures/refit_revSuppfigures2.pdf",
  plot = finals2,
  width = 12,
  height = 10,
  units = "in",
  device = cairo_pdf
)


ggsave(
  filename = "outputs/figures/refit_revSuppfigures2.tiff",
  plot = finals2,  
  width = 12,
  height = 10,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)


ggsave(
  filename = "outputs/figures/refit_Suppfigures3.tiff",
  plot = finals3,  
  width = 16,
  height = 13,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

ggsave(
  filename = "outputs/figures/refit_Suppfigures3.pdf",
  plot = finals3,
  width = 16,
  height = 13,
  units = "in",
  device = cairo_pdf
)

#ESS estimates#
# Refit version: the refitted draws are thinned to every 20th sample, so R-hat
# and ESS are not recomputed from them. They were computed on the full chains
# during extraction and are read from the scalar summaries instead.
params_to_plot <- c(
  "log_mult","transit_time_mean","transit_time_cv","tau_ww",
  "beta","kappa","report_frac","phi",
  "delta_inv","theta_invall","omega_invall","Vea","Veb","m"
)

scalar_cas <- as.data.frame(read_excel(file.path("outputs", "summaries", "refitted_case_scalar_summary.xlsx")))
scalar_ww  <- as.data.frame(read_excel(file.path("outputs", "summaries", "refitted_ww_scalar_summary.xlsx")))
scalar_com <- as.data.frame(read_excel(file.path("outputs", "summaries", "refitted_combined_scalar_summary.xlsx")))

# Check the Excel files have come through intact before using them
cat("\n--- combined model, psi rows from the scalar summary ---\n")
print(scalar_com[scalar_com$parameter %in% c("log_mult", "mult"), ], row.names = FALSE)
cat("\n--- combined model, Rhat column from the scalar summary ---\n")
print(scalar_com[, c("parameter", "Rhat")], row.names = FALSE)

# wide table (parameters as rows, models as columns)
wide_rhat <- data.frame(
  Parameter   = params_to_plot,
  `Case-only` = scalar_cas$Rhat[match(params_to_plot, scalar_cas$parameter)],
  `WW-only`   = scalar_ww$Rhat[match(params_to_plot, scalar_ww$parameter)],
  Combined    = scalar_com$Rhat[match(params_to_plot, scalar_com$parameter)],
  check.names = FALSE
)

wide_ess <- data.frame(
  Parameter   = params_to_plot,
  `Case-only` = scalar_cas$ESS[match(params_to_plot, scalar_cas$parameter)],
  `WW-only`   = scalar_ww$ESS[match(params_to_plot, scalar_ww$parameter)],
  Combined    = scalar_com$ESS[match(params_to_plot, scalar_com$parameter)],
  check.names = FALSE
)

cat("\n--- R-hat and ESS, all three refitted models ---\n")
print(wide_rhat, row.names = FALSE)
print(wide_ess, row.names = FALSE)
cat(sprintf("maximum R-hat: %.4f   minimum ESS: %.0f\n",
            max(as.matrix(wide_rhat[, -1]), na.rm = TRUE),
            min(as.matrix(wide_ess[, -1]), na.rm = TRUE)))

#  write to CSV
write.csv(wide_ess, "outputs/tables/refit_wide_ess.csv", row.names = FALSE)


# ==================================================
# Under-reporting factor (UF) estimation
# ==================================================

# Convert fitted objects to matrices
# Refit version: mcmc_matrixallcom and mcmc_matrixallcas already hold the refitted
# draws from the top of the script.
mcmc_matrixallww  <- mcmc_matrixallWW

post <- mcmc_matrixallcom   # Posterior samples (burn-in refers to time, not iterations)


# ==================================================
# 1. Symptomatic infections (model-based UF)
# ==================================================

# Extract symptomatic incidence columns
sympt_cols <- grep("^sympt_new", colnames(post))
sympt_mat  <- post[, sympt_cols]

# Remove first 30 timepoints (burn-in period)
sympt_mat_clean <- sympt_mat[, -c(1:30)]

# Total symptomatic infections per posterior draw
total_sympt <- rowSums(sympt_mat_clean, na.rm = TRUE)

# Load observed case data
case_dat      <- read_excel(file.path("data", "raw", "case_data_V2.xlsx"), sheet = "cases")
cases_obsb    <- case_dat$total_cases
cum_reported_cases <- sum(cases_obsb)

# Under-reporting factor (symptomatic-based)
UF_sympt_data <- total_sympt / cum_reported_cases

UF_summarysymp <- c(
  median   = median(UF_sympt_data),
  lower_95 = quantile(UF_sympt_data, 0.025),
  upper_95 = quantile(UF_sympt_data, 0.975)
)

UF_summarysymp


# ==================================================
# 2. Symptomatic UF — Case-only model
# ==================================================

sympt_cols <- grep("^sympt_new", colnames(mcmc_matrixallcas))
sympt_matcas <- mcmc_matrixallcas[, sympt_cols]
sympt_mat_cleancas <- sympt_matcas[, -c(1:30)]

total_symptcas <- rowSums(sympt_mat_cleancas, na.rm = TRUE)
UF_sympt_cas   <- total_symptcas / cum_reported_cases

UF_summarysympcas <- c(
  median   = median(UF_sympt_cas),
  lower_95 = quantile(UF_sympt_cas, 0.025),
  upper_95 = quantile(UF_sympt_cas, 0.975)
)

UF_summarysympcas

# ==================================================
# 3. Symptomatic UF — WW-only model
# ==================================================

sympt_cols <- grep("^sympt_new", colnames(mcmc_matrixallww))
sympt_matww <- mcmc_matrixallww[, sympt_cols]
sympt_mat_cleanww <- sympt_matww[, -c(1:30)]

total_symptww <- rowSums(sympt_mat_cleanww, na.rm = TRUE)
UF_sympt_ww   <- total_symptww / cum_reported_cases

UF_summarysympww <- c(
  median   = median(UF_sympt_ww),
  lower_95 = quantile(UF_sympt_ww, 0.025),
  upper_95 = quantile(UF_sympt_ww, 0.975)
)

UF_summarysympww


# ==================================================
# 4.Wastewater implied infections combined model
# ==================================================

R_cols <- grep("^R_total\\[", colnames(post))
R_mat  <- mcmc_matrixallww[, R_cols]

# Remove warm-up period
R_mat_clean <- R_mat[, -c(1:30)]

# Final epidemic size (last timepoint)
R_total_clean <- R_mat_clean[, ncol(R_mat_clean)]

UF_R <- R_total_clean / cum_reported_cases

UF_R_summary <- c(
  median   = median(UF_R),
  lower_95 = quantile(UF_R, 0.025),
  upper_95 = quantile(UF_R, 0.975)
)

UF_R_summary


