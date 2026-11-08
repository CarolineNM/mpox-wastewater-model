# ==================================================
# Script: extract_refit_summary.R
# Purpose: Posterior summaries, convergence diagnostics and thinned draws
#          from a fitted model
# Project: Mpox_clean
# Updated: 2026-11-07
# ==================================================

rm(list = ls())

# ---- Packages ----
library(coda)

options(scipen = 999)


# ---- Fit to extract ----
# combined_fit_refit.RData / combined_fit_base / combined
# ww_fit_refit.RData       / ww_fit            / ww
# case_fit_refit.RData     / case_fit          / case

fit_filename <- "combined_fit_refit.RData"
fit_object   <- "combined_fit_base"
prefix       <- "combined"

keep_every <- 20


# ---- Parameters to summarise ----
scalars <- c(
  "beta", "kappa", "phi", "report_frac", "Vea", "Veb", "m",
  "delta_inv", "theta_invall", "omega_invall",
  "log_mult", "mult", "tau_ww",
  "transit_time_mean", "transit_time_cv",
  "E0", "P0", "A10", "I10"
)

series <- c(
  "total_Cuminc", "Cuminc", "R_total", "sympt_new", "total_new_cases",
  "cases_pred", "mu_nb", "active_infected", "total_lambda",
  "P_total", "A_total", "I_total",
  "shed_P", "shed_A", "shed_I", "daily_shedding",
  "log10_conc", "ww_pred"
)


# ---- Load fitted model ----
load(file.path("outputs", "fits", fit_filename))
post <- get(fit_object)

draws   <- as.matrix(post)
present <- scalars[scalars %in% colnames(draws)]


# ---- Convergence diagnostics for Table S1 ----
rhat <- rep(NA_real_, length(present))
ess  <- rep(NA_real_, length(present))

for (i in seq_along(present)) {
  sub <- post[, present[i], drop = FALSE]
  rhat[i] <- gelman.diag(sub, autoburnin = FALSE)$psrf[1, 1]
  ess[i]  <- effectiveSize(sub)[[1]]
}


# ---- Scalar summary ----
summ <- data.frame(
  parameter = present,
  median    = apply(draws[, present], 2, median),
  lower95   = apply(draws[, present], 2, quantile, 0.025),
  upper95   = apply(draws[, present], 2, quantile, 0.975),
  mean      = apply(draws[, present], 2, mean),
  sd        = apply(draws[, present], 2, sd),
  Rhat      = round(rhat, 4),
  ESS       = round(ess, 0),
  row.names = NULL
)

write.csv(summ, file.path("outputs", "summaries",
                          paste0(prefix, "_scalar_summary.csv")),
          row.names = FALSE)

summ


# ---- Position of psi within its prior range ----
# reported in the limitations: distance from the nearer truncation bound
psi <- draws[, "mult"]

cat("psi median      :", signif(median(psi), 3), "\n")
cat("psi log10 median:", round(median(log10(psi)), 2), "\n")
cat("log10 units above lower bound (1e-8):",
    round(median(log10(psi)) - log10(1e-8), 2), "\n")


# ---- Thinned draws for the figure and validation scripts ----
keep_cols <- c(present, unlist(lapply(series, function(s)
  grep(paste0("^", s, "(\\[|$)"), colnames(draws), value = TRUE))))
keep_cols <- unique(keep_cols[keep_cols %in% colnames(draws)])

rows  <- seq(1, nrow(draws), by = keep_every)
small <- draws[rows, keep_cols]

saveRDS(small, file.path("outputs", "summaries",
                         paste0(prefix, "_key_draws.rds")))

dim(small)
