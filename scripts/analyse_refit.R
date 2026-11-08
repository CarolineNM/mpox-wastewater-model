###############################################################################
# analyse_refit.R
#
# Purpose: burden estimates and manuscript figures from the refitted models
# Project: Mpox_clean
#
# Inputs : outputs/summaries/{case,ww,combined}_key_draws.rds
#          data/raw/case_data_V2.xlsx
#
# Structure follows generate_manuscript_figures.R
#
# Run from the project root.
###############################################################################

library(tidyverse)
library(readxl)
library(patchwork)
library(ggridges)
library(truncnorm)

options(scipen = 999)
set.seed(123)

summ_dir <- file.path("outputs", "summaries")
fig_dir  <- file.path("outputs", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

burn_in_timesteps <- 30

# Axis limits. coord_cartesian clips rather than drops, so bands wider than
# these are still drawn. The case-only viral load band runs below zero and
# above 8, so the viral load range is widened to show it. The WW-only and
# combined viral load panels keep the narrower range.
CASE_YLIM <- c(0, 11)
CASE_YLIM_PPD <- c(0, 20)
WW_YLIM <- c(-3, 10)
WW_YLIM_WWCOM <- c(0, 8)


# ---- Load observed data ----
data_file <- file.path("data", "raw", "case_data_V2.xlsx")
case_dat  <- read_excel(data_file, sheet = "cases")
ww_dat    <- read_excel(data_file, sheet = "dailyWW")

cases_obsb <- case_dat$total_cases
ww_obs     <- as.numeric(unlist(ww_dat$log10_cp_per_person_per_day))

T_caseobs <- length(cases_obsb)
T_WWobs   <- length(ww_obs)


# ---- Load posterior draws ----
# Refitted models, used for the figures and for the new numbers
mcmc_matrixallcas <- readRDS(file.path(summ_dir, "refitted_case_key_draws.rds"))
mcmc_matrixallWW  <- readRDS(file.path(summ_dir, "refitted_ww_key_draws.rds"))
mcmc_matrixallcom <- readRDS(file.path(summ_dir, "refitted_combined_key_draws.rds"))

# Published models, used only for the before and after comparison
pub_cas <- readRDS(file.path(summ_dir, "published_case_key_draws.rds"))
pub_ww  <- readRDS(file.path(summ_dir, "published_ww_key_draws.rds"))
pub_com <- readRDS(file.path(summ_dir, "published_combined_key_draws.rds"))


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


# ==================================================
# Column names for the series used below
# ==================================================

obs_window <- (burn_in_timesteps + 1):(burn_in_timesteps + T_caseobs)

cas_cols  <- paste0("total_new_cases[", obs_window, "]")
cum_cols  <- paste0("total_Cuminc[", obs_window, "]")
symp_cols <- paste0("sympt_new[", obs_window, "]")
conc_cols <- paste0("log10_conc[", 1:T_WWobs, "]")
wwpd_cols <- paste0("ww_pred[", 1:T_WWobs, "]")


# ==================================================
# Posterior predictive for reported cases
# ==================================================
# cases_pred is monitored directly by the model, indexed from 1 to T_caseobs.

cpred_cols <- paste0("cases_pred[", 1:T_caseobs, "]")

casepred_cas <- mcmc_matrixallcas[, cpred_cols]
casepred_ww  <- mcmc_matrixallWW[,  cpred_cols]
casepred_com <- mcmc_matrixallcom[, cpred_cols]


# ==================================================
# Extract model-specific summaries
# ==================================================

# fitted mean trajectories
total_cas_cas_summary <- summary_median_CI(mcmc_matrixallcas[, cas_cols])
total_ww_cas_summary  <- summary_median_CI(mcmc_matrixallWW[,  cas_cols])
total_com_cas_summary <- summary_median_CI(mcmc_matrixallcom[, cas_cols])

total_cas_ww_summary  <- summary_median_CI(mcmc_matrixallcas[, conc_cols])
total_ww_ww_summary   <- summary_median_CI(mcmc_matrixallWW[,  conc_cols])
total_com_ww_summary  <- summary_median_CI(mcmc_matrixallcom[, conc_cols])

# posterior predictive
pred_cas_cas_summary <- summary_median_CI(casepred_cas)
pred_ww_cas_summary  <- summary_median_CI(casepred_ww)
pred_com_cas_summary <- summary_median_CI(casepred_com)

pred_cas_ww_summary  <- summary_median_CI(mcmc_matrixallcas[, wwpd_cols])
pred_ww_ww_summary   <- summary_median_CI(mcmc_matrixallWW[,  wwpd_cols])
pred_com_ww_summary  <- summary_median_CI(mcmc_matrixallcom[, wwpd_cols])


# ==================================================
# Plot data frames
# ==================================================

# ---- fitted means, cases ----
plot_datcas <- data.frame(
  time       = 1:T_caseobs,
  Date       = case_dat$Date,
  observed   = cases_obsb,
  median_fit = total_cas_cas_summary[, "median"],
  lower_ci   = total_cas_cas_summary[, "lower_95_CI"],
  upper_ci   = total_cas_cas_summary[, "upper_95_CI"]
)

plot_datwwcas <- data.frame(
  time       = 1:T_caseobs,
  Date       = case_dat$Date,
  observed   = cases_obsb,
  median_fit = total_ww_cas_summary[, "median"],
  lower_ci   = total_ww_cas_summary[, "lower_95_CI"],
  upper_ci   = total_ww_cas_summary[, "upper_95_CI"]
)

plot_datcomcas <- data.frame(
  time       = 1:T_caseobs,
  Date       = case_dat$Date,
  observed   = cases_obsb,
  median_fit = total_com_cas_summary[, "median"],
  lower_ci   = total_com_cas_summary[, "lower_95_CI"],
  upper_ci   = total_com_cas_summary[, "upper_95_CI"]
)

# ---- fitted means, viral load ----
plot_datcasww <- data.frame(
  time       = 1:T_WWobs,
  Date       = ww_dat$date,
  observed   = ww_obs,
  median_fit = total_cas_ww_summary[, "median"],
  lower_ci   = total_cas_ww_summary[, "lower_95_CI"],
  upper_ci   = total_cas_ww_summary[, "upper_95_CI"]
)

plot_datww <- data.frame(
  time       = 1:T_WWobs,
  Date       = ww_dat$date,
  observed   = ww_obs,
  median_fit = total_ww_ww_summary[, "median"],
  lower_ci   = total_ww_ww_summary[, "lower_95_CI"],
  upper_ci   = total_ww_ww_summary[, "upper_95_CI"]
)

plot_datcomww <- data.frame(
  time       = 1:T_WWobs,
  Date       = ww_dat$date,
  observed   = ww_obs,
  median_fit = total_com_ww_summary[, "median"],
  lower_ci   = total_com_ww_summary[, "lower_95_CI"],
  upper_ci   = total_com_ww_summary[, "upper_95_CI"]
)

# ---- posterior predictive, cases ----
ppd_datcas <- data.frame(
  time       = 1:T_caseobs,
  Date       = case_dat$Date,
  observed   = cases_obsb,
  median_fit = pred_cas_cas_summary[, "median"],
  lower_ci   = pred_cas_cas_summary[, "lower_95_CI"],
  upper_ci   = pred_cas_cas_summary[, "upper_95_CI"]
)

ppd_datwwcas <- data.frame(
  time       = 1:T_caseobs,
  Date       = case_dat$Date,
  observed   = cases_obsb,
  median_fit = pred_ww_cas_summary[, "median"],
  lower_ci   = pred_ww_cas_summary[, "lower_95_CI"],
  upper_ci   = pred_ww_cas_summary[, "upper_95_CI"]
)

ppd_datcomcas <- data.frame(
  time       = 1:T_caseobs,
  Date       = case_dat$Date,
  observed   = cases_obsb,
  median_fit = pred_com_cas_summary[, "median"],
  lower_ci   = pred_com_cas_summary[, "lower_95_CI"],
  upper_ci   = pred_com_cas_summary[, "upper_95_CI"]
)

# ---- posterior predictive, viral load ----
ppd_datcasww <- data.frame(
  time       = 1:T_WWobs,
  Date       = ww_dat$date,
  observed   = ww_obs,
  median_fit = pred_cas_ww_summary[, "median"],
  lower_ci   = pred_cas_ww_summary[, "lower_95_CI"],
  upper_ci   = pred_cas_ww_summary[, "upper_95_CI"]
)

ppd_datww <- data.frame(
  time       = 1:T_WWobs,
  Date       = ww_dat$date,
  observed   = ww_obs,
  median_fit = pred_ww_ww_summary[, "median"],
  lower_ci   = pred_ww_ww_summary[, "lower_95_CI"],
  upper_ci   = pred_ww_ww_summary[, "upper_95_CI"]
)

ppd_datcomww <- data.frame(
  time       = 1:T_WWobs,
  Date       = ww_dat$date,
  observed   = ww_obs,
  median_fit = pred_com_ww_summary[, "median"],
  lower_ci   = pred_com_ww_summary[, "lower_95_CI"],
  upper_ci   = pred_com_ww_summary[, "upper_95_CI"]
)


# ==================================================
# How far do the ribbons extend?
# ==================================================
# coord_cartesian clips rather than drops, so bands wider than the axis are
# still drawn. This prints their true extent so the axis ranges can be chosen
# deliberately rather than by accident.

ribbon_range <- function(d, label) {
  cat(sprintf("  %-34s lower %8.2f to %8.2f   upper %8.2f to %8.2f\n",
              label, min(d$lower_ci), max(d$lower_ci),
              min(d$upper_ci), max(d$upper_ci)))
}

cat("\n--- fitted means ---\n")
ribbon_range(plot_datcas,    "cases, case-only")
ribbon_range(plot_datwwcas,  "cases, WW-only")
ribbon_range(plot_datcomcas, "cases, combined")
ribbon_range(plot_datcasww,  "viral load, case-only")
ribbon_range(plot_datww,     "viral load, WW-only")
ribbon_range(plot_datcomww,  "viral load, combined")

cat("\n--- posterior predictive ---\n")
ribbon_range(ppd_datcas,    "cases, case-only")
ribbon_range(ppd_datwwcas,  "cases, WW-only")
ribbon_range(ppd_datcomcas, "cases, combined")
ribbon_range(ppd_datcasww,  "viral load, case-only")
ribbon_range(ppd_datww,     "viral load, WW-only")
ribbon_range(ppd_datcomww,  "viral load, combined")
cat("\n")


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
# Case fits, fitted means (main figure)
# ==================================================

plot_casefit <- ggplot(plot_datcas, aes(x = Date)) +
  plot_geom +
  coord_cartesian(ylim = CASE_YLIM) +
  labs(x = "Date", y = "Fit vs. observed cases",
       title = "(Case-only model)") +
  custom_theme

plot_wwcasefit <- ggplot(plot_datwwcas, aes(x = Date)) +
  plot_geom +
  coord_cartesian(ylim = CASE_YLIM) +
  labs(x = "Date", y = "Fit vs. observed cases",
       title = "(WW-only model)") +
  custom_theme

plot_comcasefit <- ggplot(plot_datcomcas, aes(x = Date)) +
  plot_geom +
  coord_cartesian(ylim = CASE_YLIM) +
  labs(x = "Date", y = "Fit vs. observed cases",
       title = "(Combined model)") +
  custom_theme


# ==================================================
# Viral load fits, fitted means (main figure)
# ==================================================

plot_caswwfit <- ggplot(plot_datcasww, aes(x = Date)) +
  plot_geom +
  coord_cartesian(ylim = WW_YLIM) +
  labs(x = "Date", y = "Fit vs. observed viral load",
       title = "(Case-only model)") +
  custom_theme

plot_wwfit <- ggplot(plot_datww, aes(x = Date)) +
  plot_geom +
  coord_cartesian(ylim = WW_YLIM_WWCOM) +
  labs(x = "Date", y = "Fit vs. observed viral load",
       title = "(WW-only model)") +
  custom_theme

plot_comwwfit <- ggplot(plot_datcomww, aes(x = Date)) +
  plot_geom +
  coord_cartesian(ylim = WW_YLIM_WWCOM) +
  labs(x = "Date", y = "Fit vs. observed viral load",
       title = "(Combined model)") +
  custom_theme


# ==================================================
# Posterior predictive, cases (supplementary figure S1)
# ==================================================

ppd_casefit <- ggplot(ppd_datcas, aes(x = Date)) +
  plot_geom +
  coord_cartesian(ylim = CASE_YLIM_PPD) +
  labs(x = "Date", y = "Fit vs. observed cases",
       title = "(Case-only model)") +
  custom_theme

ppd_wwcasefit <- ggplot(ppd_datwwcas, aes(x = Date)) +
  plot_geom +
  coord_cartesian(ylim = CASE_YLIM_PPD) +
  labs(x = "Date", y = "Fit vs. observed cases",
       title = "(WW-only model)") +
  custom_theme

ppd_comcasefit <- ggplot(ppd_datcomcas, aes(x = Date)) +
  plot_geom +
  coord_cartesian(ylim = CASE_YLIM_PPD) +
  labs(x = "Date", y = "Fit vs. observed cases",
       title = "(Combined model)") +
  custom_theme


# ==================================================
# Posterior predictive, viral load (supplementary figure S1)
# ==================================================

ppd_caswwfit <- ggplot(ppd_datcasww, aes(x = Date)) +
  plot_geom +
  coord_cartesian(ylim = WW_YLIM) +
  labs(x = "Date", y = "Fit vs. observed viral load",
       title = "(Case-only model)") +
  custom_theme

ppd_wwfit <- ggplot(ppd_datww, aes(x = Date)) +
  plot_geom +
  coord_cartesian(ylim = WW_YLIM_WWCOM) +
  labs(x = "Date", y = "Fit vs. observed viral load",
       title = "(WW-only model)") +
  custom_theme

ppd_comwwfit <- ggplot(ppd_datcomww, aes(x = Date)) +
  plot_geom +
  coord_cartesian(ylim = WW_YLIM_WWCOM) +
  labs(x = "Date", y = "Fit vs. observed viral load",
       title = "(Combined model)") +
  custom_theme


# ==================================================
# Combine and save
# ==================================================

# Main figure: fitted means
top_row    <- plot_casefit + plot_wwcasefit + plot_comcasefit
bottom_row <- plot_caswwfit + plot_wwfit + plot_comwwfit

mainfigure1 <- (top_row / bottom_row) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Supplementary figure S1: posterior predictive
top_rowb    <- ppd_casefit + ppd_wwcasefit + ppd_comcasefit
bottom_rowb <- ppd_caswwfit + ppd_wwfit + ppd_comwwfit

figures1 <- (top_rowb / bottom_rowb) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = "outputs/figures/diag_mainfigure1.pdf",
  plot = mainfigure1,
  width = 12, height = 7, units = "in", device = cairo_pdf
)

ggsave(
  filename = "outputs/figures/diag_mainfigure1.tiff",
  plot = mainfigure1,
  width = 12, height = 7, dpi = 300, units = "in",
  device = "tiff", compression = "lzw"
)

ggsave(
  filename = "outputs/figures/diag_suppfigure1.pdf",
  plot = figures1,
  width = 12, height = 7, units = "in", device = cairo_pdf
)

ggsave(
  filename = "outputs/figures/diag_suppfigure1.tiff",
  plot = figures1,
  width = 12, height = 7, dpi = 300, units = "in",
  device = "tiff", compression = "lzw"
)


# ==================================================
# Supplementary figure S2: prior versus posterior
# ==================================================
# NOTE: the prior on the wastewater scaling parameter has changed. The
# definition below is the revised one. The previous definition was
#   exp(rtruncnorm(n, log(1e-9), log(1e-8), log(3e-9), sqrt(1/2.5)))
# and must not be used with the refitted posteriors.
# psi now spans eight orders of magnitude, so it is shown on the log10 scale.

prior_definitions <- list(
  kappa             = function(n) rbeta(n, 40, 2),
  report_frac       = function(n) rbeta(n, 10, 10),
  beta              = function(n) rtruncnorm(n, 0, 1, mean = 0.8, sd = 0.1),
  log10_mult        = function(n) log10(exp(rtruncnorm(n, a = log(1e-8), b = log(1e0),
                                                       mean = log(1.5e-4), sd = sqrt(1 / 0.05)))),
  transit_time_mean = function(n) rtruncnorm(n, a = 0.1, b = 10, mean = 2.5, sd = sqrt(1 / 0.25)),
  tau_ww            = function(n) rgamma(n, shape = 40, rate = 48)
)

param_labels <- c(
  kappa             = "Mixing probability, κ",
  report_frac       = "Reporting fraction",
  beta              = "Transmission probability, β",
  log10_mult        = "log10(ψ)",
  transit_time_mean = "Mean transit time (days)",
  tau_ww            = "Precision of WW likelihood, τ_ww"
)

TR_params <- c("kappa", "report_frac", "beta")
WW_params <- c("log10_mult", "transit_time_mean", "tau_ww")

n_samples <- 40000
all_param_dfs <- list()

for (param in names(param_labels)) {
  col <- if (param == "log10_mult") "mult" else param

  prior_df <- data.frame(
    value     = prior_definitions[[param]](n_samples),
    source    = "Prior",
    parameter = param
  )

  v_cas <- mcmc_matrixallcas[, col]
  v_ww  <- mcmc_matrixallWW[,  col]
  v_com <- mcmc_matrixallcom[, col]

  if (param == "log10_mult") {
    v_cas <- log10(v_cas); v_ww <- log10(v_ww); v_com <- log10(v_com)
  }

  all_param_dfs[[param]] <- bind_rows(
    prior_df,
    data.frame(value = v_cas, source = "Case-only",       parameter = param),
    data.frame(value = v_ww,  source = "Wastewater-only", parameter = param),
    data.frame(value = v_com, source = "Combined",        parameter = param)
  )
}

beta_df_all <- bind_rows(all_param_dfs)

# ggridges draws the first level at the bottom, so Prior appears at the top
beta_df_all$source <- factor(
  beta_df_all$source,
  levels = c("Combined", "Wastewater-only", "Case-only", "Prior")
)

beta_df_all$parameter_label <- factor(
  beta_df_all$parameter,
  levels = names(param_labels),
  labels = param_labels
)

grey_fills <- c("Prior"           = "grey95",
                "Case-only"       = "grey85",
                "Wastewater-only" = "grey70",
                "Combined"        = "grey45")

plota <- ggplot(beta_df_all %>% filter(parameter %in% TR_params),
                aes(x = value, y = source, fill = source)) +
  geom_density_ridges(alpha = 0.9, quantile_lines = TRUE, quantiles = 0.5,
                      scale = 0.9, rel_min_height = 0.0005,
                      colour = "grey20", linewidth = 0.3) +
  facet_wrap(~ parameter_label, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = grey_fills) +
  labs(x = "Parameter value", y = "Source",
       title = "Transmission and surveillance parameters") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5, colour = "black"),
        strip.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 11, face = "bold", colour = "black"),
        axis.text  = element_text(size = 9, colour = "black"))

plotb <- ggplot(beta_df_all %>% filter(parameter %in% WW_params),
                aes(x = value, y = source, fill = source)) +
  geom_density_ridges(alpha = 0.9, quantile_lines = TRUE, quantiles = 0.5,
                      scale = 0.9, rel_min_height = 0.0005,
                      colour = "grey20", linewidth = 0.3) +
  facet_wrap(~ parameter_label, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = grey_fills) +
  labs(x = "Parameter value", y = "Source",
       title = "Wastewater parameters") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5, colour = "black"),
        strip.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 11, face = "bold", colour = "black"),
        axis.text  = element_text(size = 9, colour = "black"))

suppfigure2 <- (plota / plotb) +
  plot_annotation(tag_levels = "A")

ggsave(
  filename = "outputs/figures/diag_suppfigure2.pdf",
  plot = suppfigure2,
  width = 10, height = 9, units = "in", device = cairo_pdf
)

ggsave(
  filename = "outputs/figures/diag_suppfigure2.tiff",
  plot = suppfigure2,
  width = 10, height = 9, dpi = 300, units = "in",
  device = "tiff", compression = "lzw"
)


# ==================================================
# Burden and underascertainment
# ==================================================
# The denominator of each factor can be read as the observed cumulative
# reported cases or as the model-predicted cumulative reported cases. Both are
# reported so the definition reproducing the published values can be identified.

# Definitions follow generate_manuscript_figures.R: symptomatic infections
# summed over the observation window, and cumulative infections taken as
# R_total at the final timepoint, each divided by the OBSERVED cumulative
# reported cases.

total_observed_cases <- sum(cases_obsb)

rtot_cols <- paste0("R_total[", obs_window, "]")
last_rtot <- rtot_cols[length(rtot_cols)]
last_cum  <- cum_cols[length(cum_cols)]
first_cum <- paste0("total_Cuminc[", burn_in_timesteps, "]")

fmt <- function(x, digits = 0) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f to %.", digits, "f)"),
          median(x), quantile(x, 0.025), quantile(x, 0.975))
}

# ---- published, case-only ----
pcas_cum  <- pub_cas[, last_cum] - pub_cas[, first_cum]
pcas_inf  <- pub_cas[, last_rtot]
pcas_symp <- rowSums(pub_cas[, symp_cols])

# ---- published, wastewater-only ----
pww_cum  <- pub_ww[, last_cum] - pub_ww[, first_cum]
pww_inf  <- pub_ww[, last_rtot]
pww_symp <- rowSums(pub_ww[, symp_cols])

# ---- published, combined ----
pcom_cum  <- pub_com[, last_cum] - pub_com[, first_cum]
pcom_inf  <- pub_com[, last_rtot]
pcom_symp <- rowSums(pub_com[, symp_cols])

# ---- refitted, case-only ----
rcas_cum  <- mcmc_matrixallcas[, last_cum] - mcmc_matrixallcas[, first_cum]
rcas_inf  <- mcmc_matrixallcas[, last_rtot]
rcas_symp <- rowSums(mcmc_matrixallcas[, symp_cols])

# ---- refitted, wastewater-only ----
rww_cum  <- mcmc_matrixallWW[, last_cum] - mcmc_matrixallWW[, first_cum]
rww_inf  <- mcmc_matrixallWW[, last_rtot]
rww_symp <- rowSums(mcmc_matrixallWW[, symp_cols])

# ---- refitted, combined ----
rcom_cum  <- mcmc_matrixallcom[, last_cum] - mcmc_matrixallcom[, first_cum]
rcom_inf  <- mcmc_matrixallcom[, last_rtot]
rcom_symp <- rowSums(mcmc_matrixallcom[, symp_cols])

burden <- data.frame(

  fit = c("Published", "Published", "Published",
          "Refitted",  "Refitted",  "Refitted"),

  model = c("Case-only", "Wastewater-only", "Combined",
            "Case-only", "Wastewater-only", "Combined"),

  cumulative_incidence = c(fmt(pcas_cum), fmt(pww_cum), fmt(pcom_cum),
                           fmt(rcas_cum), fmt(rww_cum), fmt(rcom_cum)),

  total_infections_R = c(fmt(pcas_inf), fmt(pww_inf), fmt(pcom_inf),
                         fmt(rcas_inf), fmt(rww_inf), fmt(rcom_inf)),

  symptomatic_factor = c(fmt(pcas_symp / total_observed_cases, 2),
                         fmt(pww_symp  / total_observed_cases, 2),
                         fmt(pcom_symp / total_observed_cases, 2),
                         fmt(rcas_symp / total_observed_cases, 2),
                         fmt(rww_symp  / total_observed_cases, 2),
                         fmt(rcom_symp / total_observed_cases, 2)),

  total_factor = c(fmt(pcas_inf / total_observed_cases, 2),
                   fmt(pww_inf  / total_observed_cases, 2),
                   fmt(pcom_inf / total_observed_cases, 2),
                   fmt(rcas_inf / total_observed_cases, 2),
                   fmt(rww_inf  / total_observed_cases, 2),
                   fmt(rcom_inf / total_observed_cases, 2)),

  transit_time = c(fmt(pub_cas[, "transit_time_mean"], 2),
                   fmt(pub_ww[,  "transit_time_mean"], 2),
                   fmt(pub_com[, "transit_time_mean"], 2),
                   fmt(mcmc_matrixallcas[, "transit_time_mean"], 2),
                   fmt(mcmc_matrixallWW[,  "transit_time_mean"], 2),
                   fmt(mcmc_matrixallcom[, "transit_time_mean"], 2)),

  psi_median = c(signif(median(pub_cas[, "mult"]), 4),
                 signif(median(pub_ww[,  "mult"]), 4),
                 signif(median(pub_com[, "mult"]), 4),
                 signif(median(mcmc_matrixallcas[, "mult"]), 4),
                 signif(median(mcmc_matrixallWW[,  "mult"]), 4),
                 signif(median(mcmc_matrixallcom[, "mult"]), 4))
)

print(burden, width = 300)
write.csv(burden, file.path(summ_dir, "refit_burden_summary.csv"), row.names = FALSE)

cat("\nObserved cumulative reported cases:", total_observed_cases, "\n")


# ==================================================
# Do the log_mult excursions change anything?
# ==================================================
# A small fraction of draws in the refitted combined model run up toward the
# upper bound of the psi prior. This checks whether excluding them moves any
# reported number.

keep <- mcmc_matrixallcom[, "log_mult"] <= -6

cat("\n--- log_mult excursions, refitted combined ---\n")
cat(sprintf("  draws with log_mult > -6 : %d of %d  (%.2f%%)\n",
            sum(!keep), length(keep), 100 * mean(!keep)))

cat(sprintf("  %-22s all draws: %-28s excluding excursions: %s\n",
            "cumulative incidence", fmt(rcom_cum), fmt(rcom_cum[keep])))
cat(sprintf("  %-22s all draws: %-28s excluding excursions: %s\n",
            "total infections", fmt(rcom_inf), fmt(rcom_inf[keep])))
cat(sprintf("  %-22s all draws: %-28s excluding excursions: %s\n",
            "symptomatic factor",
            fmt(rcom_symp / total_observed_cases, 2),
            fmt(rcom_symp[keep] / total_observed_cases, 2)))
cat(sprintf("  %-22s all draws: %-28s excluding excursions: %s\n",
            "total factor",
            fmt(rcom_inf / total_observed_cases, 2),
            fmt(rcom_inf[keep] / total_observed_cases, 2)))
cat(sprintf("  %-22s all draws: %-28s excluding excursions: %s\n",
            "psi",
            sprintf("%.4g", median(mcmc_matrixallcom[, "mult"])),
            sprintf("%.4g", median(mcmc_matrixallcom[keep, "mult"]))))

cat("\nFigures written to ", fig_dir, "\n", sep = "")
