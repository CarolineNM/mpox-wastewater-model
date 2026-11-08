# ==================================================
# Script: extract_extra_series.R
# Purpose: Group-level wastewater concentration and expected case series,
#          used by panel B of manuscript figure 4
# Project: Mpox_clean
# Updated: 2026-11-07
# ==================================================

rm(list = ls())

options(scipen = 999)


# ---- Fit to extract ----
# Thinning must match extract_refit_summary.R so the rows line up with
# <prefix>_key_draws.rds

fit_filename <- "combined_fit_refit.RData"
fit_object   <- "combined_fit_base"
prefix       <- "refitted_combined"

keep_every <- 20

extra <- c("log10_conc_all_g", "mu_nb_g")


# ---- Load fitted model ----
load(file.path("outputs", "fits", fit_filename))

draws <- as.matrix(get(fit_object))

keep_cols <- unlist(lapply(extra, function(s)
  grep(paste0("^", s, "\\["), colnames(draws), value = TRUE)))

rows  <- seq(1, nrow(draws), by = keep_every)
small <- draws[rows, keep_cols]


# ---- Check the rows match the key draws before saving ----
key <- readRDS(file.path("outputs", "summaries",
                         paste0(prefix, "_key_draws.rds")))

cat("key draws rows:", nrow(key), "  this object rows:", nrow(small), "\n")


saveRDS(small, file.path("outputs", "summaries",
                         paste0(prefix, "_extra_draws.rds")))

dim(small)
