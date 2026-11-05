###############################################################################
# extract_sensitivity.R
#
# Pulls everything needed for Table S2, Table S3 and the decay-sensitivity
# paragraph out of the five sensitivity fits, plus the baseline combined refit.
#
# Writes ONE small csv that can be downloaded on its own, and one thinned draws
# file per fit in case anything else is needed later.
#
# Run from the project root, after the fits have finished.
#
# The baseline fit is included deliberately. Its symptomatic and total
# underascertainment factors are already known to be 5.39 and 8.97. The script
# recomputes them here and STOPS if they do not come back, because that would
# mean this script computes them differently from analyse_refit.R and the
# sensitivity rows would not be comparable to the baseline row.
###############################################################################

fit_dir <- file.path("outputs", "fits")
out_dir <- file.path("outputs", "summaries")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

KEEP_EVERY <- 20

# cumulative reported cases over the study period, the denominator of both
# underascertainment factors
OBS_CASES <- 145

# first timestep after the 30-step initialisation
T_START <- 31
T_END   <- 199

# the fits, in the order they will appear in the output table
fit_files <- c(
  "combined_fit_refit.RData",
  "combined_decay_low_refit.RData",
  "combined_decay_high_refit.RData",
  "combined_fit_shed20_refit.RData",
  "combined_fit_shed50_refit.RData",
  "combined_fit_flatshed_refit.RData"
)

fit_labels <- c(
  "Baseline (mu = 0.143)",
  "Viral decay mu = 0.10",
  "Viral decay mu = 0.20",
  "Shedding reduced 20% (shed_mult = 0.8)",
  "Shedding reduced 50% (shed_mult = 0.5)",
  "Flat shedding profile"
)

# nodes kept in the thinned draws files
series  <- c("total_Cuminc", "active_infected", "sympt_new", "R_total")
scalars <- c("beta", "kappa", "report_frac", "m", "mult", "log_mult",
             "transit_time_mean", "transit_time_cv", "tau_ww")

# ---- results table, filled in as we go -------------------------------------

res <- data.frame(
  scenario        = character(0),
  file            = character(0),
  peak_prev       = character(0),
  cum_inc         = character(0),
  sympt_UF        = character(0),
  total_UF        = character(0),
  transit_days    = character(0),
  psi             = character(0),
  max_rhat        = numeric(0),
  min_ess         = numeric(0),
  stringsAsFactors = FALSE
)

# ---- loop over the fits ----------------------------------------------------

for (k in seq_along(fit_files)) {

  fit_file <- file.path(fit_dir, fit_files[k])

  if (!file.exists(fit_file)) {
    message("MISSING: ", fit_file, "  -- skipped")
    next
  }

  message("\n=== ", fit_labels[k], " ===")
  message("Loading ", fit_file, " ...")

  env    <- new.env(parent = emptyenv())
  loaded <- load(fit_file, envir = env)

  obj_name <- loaded[vapply(loaded,
                            function(n) inherits(env[[n]], "mcmc.list"),
                            logical(1))]
  if (length(obj_name) != 1L) {
    stop("Expected exactly one mcmc.list in ", fit_file, "; found ",
         length(obj_name), call. = FALSE)
  }
  post <- env[[obj_name]]

  draws <- do.call(rbind, lapply(post, as.matrix))
  message("  object: ", obj_name, "   draws: ", nrow(draws),
          "   nodes: ", ncol(draws))

  # ---- peak prevalence: max over the observation window, per draw ----------

  prev_cols <- paste0("active_infected[", T_START:T_END, "]")
  prev_cols <- prev_cols[prev_cols %in% colnames(draws)]
  peak_draw <- apply(draws[, prev_cols, drop = FALSE], 1, max)

  # ---- cumulative incidence: total_Cuminc at the final timestep ------------

  cum_draw <- draws[, paste0("total_Cuminc[", T_END, "]")]

  # ---- underascertainment --------------------------------------------------
  # symptomatic: incident symptomatic infections summed over the observation
  # window, divided by cumulative reported cases
  # total: all infections (R_total at the final timestep) over the same

  sympt_cols <- paste0("sympt_new[", T_START:T_END, "]")
  sympt_cols <- sympt_cols[sympt_cols %in% colnames(draws)]
  sympt_draw <- rowSums(draws[, sympt_cols, drop = FALSE]) / OBS_CASES

  total_draw <- draws[, paste0("R_total[", T_END, "]")] / OBS_CASES

  # ---- scalars -------------------------------------------------------------

  transit_draw <- draws[, "transit_time_mean"]
  psi_draw     <- draws[, "mult"]

  # ---- convergence, on the full chains ------------------------------------

  present <- scalars[scalars %in% colnames(draws)]
  rhat <- rep(NA_real_, length(present))
  ess  <- rep(NA_real_, length(present))

  if (requireNamespace("coda", quietly = TRUE) && length(post) > 1) {
    for (i in seq_along(present)) {
      sub <- coda::as.mcmc.list(lapply(post, function(ch)
        coda::as.mcmc(ch[, present[i], drop = FALSE])))
      rhat[i] <- tryCatch(coda::gelman.diag(sub, autoburnin = FALSE)$psrf[1, 1],
                          error = function(e) NA_real_)
      ess[i]  <- tryCatch(coda::effectiveSize(sub)[[1]], error = function(e) NA_real_)
    }
  }

  # mult is excluded from the convergence check: R-hat on the exponentiated
  # node is an artefact, log_mult is the sampled node
  keep_diag <- present != "mult"

  # ---- format --------------------------------------------------------------

  fmt0 <- function(x) sprintf("%.0f (%.0f-%.0f)", median(x),
                              quantile(x, 0.025), quantile(x, 0.975))
  fmt2 <- function(x) sprintf("%.2f (%.2f-%.2f)", median(x),
                              quantile(x, 0.025), quantile(x, 0.975))
  fmte <- function(x) sprintf("%.3g (%.3g-%.3g)", median(x),
                              quantile(x, 0.025), quantile(x, 0.975))

  res <- rbind(res, data.frame(
    scenario     = fit_labels[k],
    file         = fit_files[k],
    peak_prev    = fmt0(peak_draw),
    cum_inc      = fmt0(cum_draw),
    sympt_UF     = fmt2(sympt_draw),
    total_UF     = fmt2(total_draw),
    transit_days = fmt2(transit_draw),
    psi          = fmte(psi_draw),
    max_rhat     = round(max(rhat[keep_diag], na.rm = TRUE), 4),
    min_ess      = round(min(ess[keep_diag], na.rm = TRUE), 0),
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  peak prevalence   : %s\n", fmt0(peak_draw)))
  cat(sprintf("  cumulative inc    : %s\n", fmt0(cum_draw)))
  cat(sprintf("  symptomatic UF    : %s\n", fmt2(sympt_draw)))
  cat(sprintf("  total UF          : %s\n", fmt2(total_draw)))
  cat(sprintf("  transit time      : %s\n", fmt2(transit_draw)))
  cat(sprintf("  psi               : %s\n", fmte(psi_draw)))
  cat(sprintf("  max R-hat         : %.4f   min ESS: %.0f\n",
              max(rhat[keep_diag], na.rm = TRUE),
              min(ess[keep_diag], na.rm = TRUE)))

  # ---- the check, on the baseline only ------------------------------------

  if (k == 1) {
    cat("\n  --- verification against analyse_refit.R ---\n")
    cat(sprintf("  symptomatic UF  expected 5.39, got %.2f\n", median(sympt_draw)))
    cat(sprintf("  total UF        expected 8.97, got %.2f\n", median(total_draw)))
    ok <- abs(median(sympt_draw) - 5.39) < 0.05 &&
          abs(median(total_draw) - 8.97) < 0.05
    if (!ok) {
      cat("\n  The baseline does not reproduce. Do not use these numbers.\n")
      cat("  Either OBS_CASES (currently ", OBS_CASES,
          ") is not the denominator analyse_refit.R uses, or the window\n",
          "  T_START:T_END is different. Check analyse_refit.R and report\n",
          "  before continuing.\n", sep = "")
      stop("baseline verification failed", call. = FALSE)
    }
    cat("  OK: this script matches analyse_refit.R.\n")
  }

  # ---- thinned draws, in case anything else is needed later ---------------

  keep_cols <- c(present,
                 unlist(lapply(series, function(s)
                   grep(paste0("^", s, "(\\[|$)"), colnames(draws), value = TRUE))))
  keep_cols <- unique(keep_cols[keep_cols %in% colnames(draws)])

  rows  <- seq(1, nrow(draws), by = KEEP_EVERY)
  small <- draws[rows, keep_cols, drop = FALSE]

  prefix <- sub("\\.RData$", "", fit_files[k])
  saveRDS(small, file.path(out_dir, paste0("sens_", prefix, "_draws.rds")))

  rm(env, post, draws, small)
  gc()
}

# ---- write the table -------------------------------------------------------

write.csv(res, file.path(out_dir, "sensitivity_summary.csv"), row.names = FALSE)

cat("\n\n=== sensitivity_summary.csv ===\n")
print(res[, c("scenario", "peak_prev", "cum_inc", "sympt_UF", "total_UF",
              "transit_days", "psi")], row.names = FALSE)

cat("\n--- written ---\n")
cat("  ", file.path(out_dir, "sensitivity_summary.csv"),
    "   <- download this one\n", sep = "")
cat("  sens_*_draws.rds, one per fit, only needed if something else comes up\n")
