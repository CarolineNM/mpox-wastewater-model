###############################################################################
# extract_refit_summary.R
#
# Purpose: pull the quantities needed from a fitted model into small files that
#          can be moved off the tower, instead of transferring the full
#          multi-gigabyte .RData objects.

###############################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript extract_refit_summary.R <fit file> <output prefix>",
       call. = FALSE)
}

fit_file <- file.path("outputs", "fits", args[1])
prefix   <- args[2]
out_dir  <- file.path("outputs", "summaries")

# keep every KEEP_EVERY-th draw in the saved draws object; the summaries below
# are computed on all draws regardless
KEEP_EVERY <- 10

# scalar parameters
scalars <- c(
  "beta", "kappa", "phi", "report_frac", "Vea", "Veb", "m",
  "delta_inv", "theta_invall", "omega_invall",
  "log_mult", "mult", "tau_ww",
  "transit_time_mean", "transit_time_cv",
  "E0", "P0", "A10", "I10"
)

# series needed to recompute burden and the underascertainment factors, and to
# check where the wastewater fit sits
series <- c("total_Cuminc", "sympt_new", "total_new_cases",
            "active_infected", "daily_shedding", "log10_conc", "ww_pred")

# bounds of the wastewater scaling prior, for the diagnostic below
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

# one matrix with all draws from all chains
draws <- do.call(rbind, lapply(post, as.matrix))
message("  total draws: ", nrow(draws), "   monitored nodes: ", ncol(draws))

# ---- scalar summary --------------------------------------------------------

present <- scalars[scalars %in% colnames(draws)]
missing <- setdiff(scalars, present)
if (length(missing) > 0) {
  message("  note: not monitored in this fit: ", paste(missing, collapse = ", "))
}

summ <- data.frame(
  parameter = present,
  median    = apply(draws[, present, drop = FALSE], 2, median),
  lower95   = apply(draws[, present, drop = FALSE], 2, quantile, 0.025),
  upper95   = apply(draws[, present, drop = FALSE], 2, quantile, 0.975),
  mean      = apply(draws[, present, drop = FALSE], 2, mean),
  sd        = apply(draws[, present, drop = FALSE], 2, sd),
  row.names = NULL
)

write.csv(summ, file.path(out_dir, paste0(prefix, "_scalar_summary.csv")),
          row.names = FALSE)

# ---- where does the wastewater scaling parameter sit? ----------------------
# A posterior concentrated against either bound indicates the scale is being
# set by the prior rather than by the data. This is the check that was not
# available before.

if ("mult" %in% colnames(draws)) {
  m <- draws[, "mult"]
  lo <- log10(MULT_LOWER); hi <- log10(MULT_UPPER)
  pos <- (log10(m) - lo) / (hi - lo)   # 0 at lower bound, 1 at upper
  cat("\n--- wastewater scaling parameter ---\n")
  cat(sprintf("  prior range      : %g to %g  (%.0f log10 units)\n",
              MULT_LOWER, MULT_UPPER, hi - lo))
  cat(sprintf("  posterior median : %.4g   (log10 = %.3f)\n",
              median(m), median(log10(m))))
  cat(sprintf("  95%% CrI          : %.4g to %.4g\n",
              quantile(m, 0.025), quantile(m, 0.975)))
  cat(sprintf("  position in range: median at %.1f%% of the way up\n",
              100 * median(pos)))
  cat(sprintf("  draws in lowest 10%% of range : %.1f%%\n",
              100 * mean(pos <= 0.10)))
  cat(sprintf("  draws in highest 10%% of range: %.1f%%\n",
              100 * mean(pos >= 0.90)))
  if (mean(pos <= 0.10) > 0.25 || mean(pos >= 0.90) > 0.25) {
    cat("  WARNING: posterior is concentrated against a bound.\n")
    cat("           Do not report this fit before investigating.\n")
  } else {
    cat("  OK: posterior lies in the interior of the prior range.\n")
  }
}

# ---- final cumulative incidence --------------------------------------------

cum_cols <- grep("^total_Cuminc\\[", colnames(draws), value = TRUE)
if (length(cum_cols) > 0) {
  idx <- as.integer(sub(".*\\[(\\d+)\\]", "\\1", cum_cols))
  last <- cum_cols[which.max(idx)]
  cum <- draws[, last]
  cat("\n--- cumulative incidence at the final timestep ---\n")
  cat(sprintf("  median %.0f   95%% CrI %.0f to %.0f\n",
              median(cum), quantile(cum, 0.025), quantile(cum, 0.975)))
}

# ---- save thinned draws ----------------------------------------------------

keep_cols <- c(present,
               unlist(lapply(series, function(s)
                 grep(paste0("^", s, "(\\[|$)"), colnames(draws), value = TRUE))))
keep_cols <- unique(keep_cols[keep_cols %in% colnames(draws)])

rows <- seq(1, nrow(draws), by = KEEP_EVERY)
small <- draws[rows, keep_cols, drop = FALSE]

saveRDS(small, file.path(out_dir, paste0(prefix, "_key_draws.rds")))

cat("\n--- written ---\n")
cat("  ", file.path(out_dir, paste0(prefix, "_scalar_summary.csv")), "\n", sep = "")
cat("  ", file.path(out_dir, paste0(prefix, "_key_draws.rds")),
    sprintf("  (%d draws x %d columns)", nrow(small), ncol(small)), "\n", sep = "")
cat(sprintf("     approximately %.1f MB\n",
            file.size(file.path(out_dir, paste0(prefix, "_key_draws.rds"))) / 1e6))
