###############################################################################
# extract_refit_summary.R
#
# Purpose: pull everything needed for the manuscript and supplement out of a
#          fitted model, into files small enough to move off the tower.
#
# This is a SEPARATE script. Do not paste it into the fitting scripts. Run it
# after a fit has finished and been saved.
#
# It writes two files per fit:
#
#   <prefix>_scalar_summary.csv   a few kB - posterior median, 95% CrI,
#                                            R-hat and effective sample size
#                                            for every scalar parameter
#   <prefix>_key_draws.rds        ~100 MB  - thinned draws of the scalars and
#                                            of every series the figure and
#                                            validation scripts use
#
# To use: set the two values below, then run the whole script from the project
# root. Repeat for each fit.
#
###############################################################################

# ---- which fit to extract --------------------------------------------------
# Edit these two lines and run the script again for each fit.
#
#   refitted:   "combined_fit_refit.RData"  prefix "combined"
#               "ww_fit_refit.RData"        prefix "ww"
#               "case_fit_refit.RData"      prefix "casefit"
#
#   published:  "combined_fit.RData"        prefix "published_combined"
#               "ww_fit.RData"              prefix "published_ww"
#               "case_fit.RData"            prefix "published_case"

fit_filename <- "combined_fit_refit.RData"
prefix       <- "combined"

# ----------------------------------------------------------------------------

fit_file <- file.path("outputs", "fits", fit_filename)
out_dir  <- file.path("outputs", "summaries")

# keep every KEEP_EVERY-th draw in the saved draws object. The summaries,
# R-hat and effective sample sizes below are computed on ALL draws.
KEEP_EVERY <- 20

# scalar parameters
scalars <- c(
  "beta", "kappa", "phi", "report_frac", "Vea", "Veb", "m",
  "delta_inv", "theta_invall", "omega_invall",
  "log_mult", "mult", "tau_ww",
  "transit_time_mean", "transit_time_cv",
  "E0", "P0", "A10", "I10"
)

# every series used by generate_manuscript_figures.R, run_counterfactuals.R
# and SeroValidation.R
series <- c(
  "total_Cuminc",     # cumulative incidence, all groups
  "Cuminc",           # cumulative incidence by activity group (figure 6)
  "R_total",          # cumulative infections, used for the total UF
  "sympt_new",        # symptomatic incidence, used for the symptomatic UF
  "total_new_cases",  # predicted reported cases
  "cases_pred",       # posterior predictive, cases
  "mu_nb",            # negative binomial mean
  "active_infected",  # prevalence
  "total_lambda",     # force of infection
  "P_total", "A_total", "I_total",   # compartment totals by group
  "shed_P", "shed_A", "shed_I",      # shedding by compartment and group
  "daily_shedding",   # total daily shedding
  "log10_conc",       # fitted viral load at sampling days
  "ww_pred"           # posterior predictive, viral load
)

# bounds of the wastewater scaling prior, for the diagnostic below.
# Set to the prior used in the fit being extracted.
MULT_LOWER <- 1e-8
MULT_UPPER <- 1e0

# ---- load ------------------------------------------------------------------

if (!file.exists(fit_file)) stop("Cannot find ", fit_file, call. = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("Loading ", fit_file, " ...")
env <- new.env(parent = emptyenv())
loaded <- load(fit_file, envir = env)

obj_name <- loaded[vapply(loaded,
                          function(n) inherits(env[[n]], "mcmc.list"),
                          logical(1))]
if (length(obj_name) != 1L) {
  stop("Expected exactly one mcmc.list in ", fit_file, "; found ",
       length(obj_name), call. = FALSE)
}
post <- env[[obj_name]]
message("  object: ", obj_name, "   chains: ", length(post),
        "   draws per chain: ", nrow(post[[1]]))

draws <- do.call(rbind, lapply(post, as.matrix))
message("  total draws: ", nrow(draws), "   monitored nodes: ", ncol(draws))

present <- scalars[scalars %in% colnames(draws)]
missing <- setdiff(scalars, present)
if (length(missing) > 0) {
  message("  note: not monitored in this fit: ", paste(missing, collapse = ", "))
}

# ---- convergence diagnostics for table S1 ----------------------------------
# Computed here, on the full chains, so the big object is not needed again.

rhat <- rep(NA_real_, length(present))
ess  <- rep(NA_real_, length(present))

if (requireNamespace("coda", quietly = TRUE) && length(post) > 1) {
  for (i in seq_along(present)) {
    p <- present[i]
    sub <- coda::as.mcmc.list(lapply(post, function(ch) coda::as.mcmc(ch[, p, drop = FALSE])))
    rhat[i] <- tryCatch(coda::gelman.diag(sub, autoburnin = FALSE)$psrf[1, 1],
                        error = function(e) NA_real_)
    ess[i]  <- tryCatch(coda::effectiveSize(sub)[[1]], error = function(e) NA_real_)
  }
} else {
  message("  note: coda unavailable or single chain; R-hat and ESS not computed")
}

# ---- scalar summary --------------------------------------------------------

summ <- data.frame(
  parameter = present,
  median    = apply(draws[, present, drop = FALSE], 2, median),
  lower95   = apply(draws[, present, drop = FALSE], 2, quantile, 0.025),
  upper95   = apply(draws[, present, drop = FALSE], 2, quantile, 0.975),
  mean      = apply(draws[, present, drop = FALSE], 2, mean),
  sd        = apply(draws[, present, drop = FALSE], 2, sd),
  Rhat      = round(rhat, 4),
  ESS       = round(ess, 0),
  row.names = NULL
)

write.csv(summ, file.path(out_dir, paste0(prefix, "_scalar_summary.csv")),
          row.names = FALSE)

cat("\n--- convergence ---\n")
if (all(is.na(summ$Rhat))) {
  cat("  not computed\n")
} else {
  cat(sprintf("  maximum R-hat : %.4f\n", max(summ$Rhat, na.rm = TRUE)))
  cat(sprintf("  minimum ESS   : %.0f\n", min(summ$ESS, na.rm = TRUE)))
  bad <- summ$parameter[!is.na(summ$Rhat) & summ$Rhat >= 1.01]
  if (length(bad) > 0) {
    cat("  WARNING: R-hat at or above 1.01 for: ", paste(bad, collapse = ", "), "\n")
  } else {
    cat("  all R-hat below 1.01\n")
  }
}

# ---- where does the wastewater scaling parameter sit? ----------------------

if ("mult" %in% colnames(draws)) {
  m <- draws[, "mult"]
  lo <- log10(MULT_LOWER); hi <- log10(MULT_UPPER)
  pos <- (log10(m) - lo) / (hi - lo)
  cat("\n--- wastewater scaling parameter ---\n")
  cat(sprintf("  prior range      : %g to %g  (%.0f log10 units)\n",
              MULT_LOWER, MULT_UPPER, hi - lo))
  cat(sprintf("  posterior median : %.4g   (log10 = %.3f)\n",
              median(m), median(log10(m))))
  cat(sprintf("  95%% CrI          : %.4g to %.4g\n",
              quantile(m, 0.025), quantile(m, 0.975)))
  cat(sprintf("  position in range: median at %.1f%% of the way up\n",
              100 * median(pos)))
  cat(sprintf("  draws in lowest 10%% of range : %.1f%%\n", 100 * mean(pos <= 0.10)))
  cat(sprintf("  draws in highest 10%% of range: %.1f%%\n", 100 * mean(pos >= 0.90)))
  if (mean(pos <= 0.10) > 0.25 || mean(pos >= 0.90) > 0.25) {
    cat("  WARNING: posterior is concentrated against a bound.\n")
  } else {
    cat("  OK: posterior lies in the interior of the prior range.\n")
  }
}

# ---- save thinned draws ----------------------------------------------------

keep_cols <- c(present,
               unlist(lapply(series, function(s)
                 grep(paste0("^", s, "(\\[|$)"), colnames(draws), value = TRUE))))
keep_cols <- unique(keep_cols[keep_cols %in% colnames(draws)])

found   <- series[vapply(series, function(s)
  any(grepl(paste0("^", s, "(\\[|$)"), colnames(draws))), logical(1))]
notfound <- setdiff(series, found)
if (length(notfound) > 0) {
  message("  note: series not monitored in this fit: ",
          paste(notfound, collapse = ", "))
}

rows  <- seq(1, nrow(draws), by = KEEP_EVERY)
small <- draws[rows, keep_cols, drop = FALSE]

saveRDS(small, file.path(out_dir, paste0(prefix, "_key_draws.rds")))

cat("\n--- written ---\n")
cat("  ", file.path(out_dir, paste0(prefix, "_scalar_summary.csv")), "\n", sep = "")
cat("  ", file.path(out_dir, paste0(prefix, "_key_draws.rds")),
    sprintf("  (%d draws x %d columns)", nrow(small), ncol(small)), "\n", sep = "")
cat(sprintf("     approximately %.0f MB\n",
            file.size(file.path(out_dir, paste0(prefix, "_key_draws.rds"))) / 1e6))
