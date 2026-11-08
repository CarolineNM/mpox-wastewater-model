# ==================================================
# Script: extract_sensitivity.R
# Purpose: Burden estimates across the decay and shedding sensitivity fits,
#          for Table S2 and Table S3
# Project: Mpox_clean
# Updated: 2026-11-07
# ==================================================

rm(list = ls())

# ---- Packages ----
library(coda)

options(scipen = 999)


# ---- Settings ----
obs_cases  <- 145          # cumulative reported cases, denominator of both UFs
t_start    <- 31           # first timestep after the 30-step initialisation
t_end      <- 199
keep_every <- 20


# ---- Fits ----
fit_files <- c(
  "combined_fit_refit.RData",
  "combined_decay_low_refit.RData",
  "combined_decay_high_refit.RData",
  "combined_fit_shed20_refit.RData",
  "combined_fit_shed50_refit.RData",
  "combined_fit_flatshed_refit.RData"
)

fit_objects <- c(
  "combined_fit_base",
  "combined_decay_low",
  "combined_decay_high",
  "combined_fit_shed20",
  "combined_fit_shed50",
  "combined_fit_base"
)

fit_labels <- c(
  "Baseline (mu = 0.143)",
  "Viral decay mu = 0.10",
  "Viral decay mu = 0.20",
  "Shedding scaled to 80% of baseline",
  "Shedding scaled to 50% of baseline",
  "Flat shedding profile"
)

scalars <- c("beta", "kappa", "report_frac", "m", "mult", "log_mult",
             "transit_time_mean", "transit_time_cv", "tau_ww")

series <- c("total_Cuminc", "active_infected", "sympt_new", "R_total")


# ---- Format posterior median and 95% CrI ----
fmt <- function(x, digits = 0) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
          median(x), quantile(x, 0.025), quantile(x, 0.975))
}


# ---- Loop over the fits ----
res <- NULL

for (k in seq_along(fit_files)) {

  load(file.path("outputs", "fits", fit_files[k]))
  post  <- get(fit_objects[k])
  draws <- as.matrix(post)

  prev_cols  <- paste0("active_infected[", t_start:t_end, "]")
  sympt_cols <- paste0("sympt_new[", t_start:t_end, "]")

  peak_draw    <- apply(draws[, prev_cols], 1, max)
  cum_draw     <- draws[, paste0("total_Cuminc[", t_end, "]")]
  sympt_draw   <- rowSums(draws[, sympt_cols]) / obs_cases
  total_draw   <- draws[, paste0("R_total[", t_end, "]")] / obs_cases
  transit_draw <- draws[, "transit_time_mean"]
  psi_draw     <- draws[, "mult"]

  # R-hat on log_mult, not mult: the diagnostic on the exponentiated node is
  # a transform artefact
  present <- scalars[scalars %in% colnames(draws)]
  present <- present[present != "mult"]

  rhat <- rep(NA_real_, length(present))
  ess  <- rep(NA_real_, length(present))

  for (i in seq_along(present)) {
    sub <- post[, present[i], drop = FALSE]
    rhat[i] <- gelman.diag(sub, autoburnin = FALSE)$psrf[1, 1]
    ess[i]  <- effectiveSize(sub)[[1]]
  }

  res <- rbind(res, data.frame(
    scenario     = fit_labels[k],
    peak_prev    = fmt(peak_draw),
    cum_inc      = fmt(cum_draw),
    sympt_UF     = fmt(sympt_draw, 2),
    total_UF     = fmt(total_draw, 2),
    transit_days = fmt(transit_draw, 2),
    psi          = signif(median(psi_draw), 2),
    max_rhat     = round(max(rhat), 4),
    min_ess      = round(min(ess), 0)
  ))

  # thinned draws, kept in case anything else is needed from these fits
  keep_cols <- c(present, unlist(lapply(series, function(s)
    grep(paste0("^", s, "(\\[|$)"), colnames(draws), value = TRUE))))
  keep_cols <- unique(keep_cols[keep_cols %in% colnames(draws)])

  rows  <- seq(1, nrow(draws), by = keep_every)
  small <- draws[rows, keep_cols]

  saveRDS(small, file.path("outputs", "summaries",
                           paste0("sens_", sub("\\.RData$", "", fit_files[k]),
                                  "_draws.rds")))

  rm(post, draws, small)
  gc()
}


# ---- Write the summary table ----
write.csv(res, file.path("outputs", "summaries", "sensitivity_summary.csv"),
          row.names = FALSE)

res
