###############################################################################
# extract_extra_series.R
#
# Pulls the two group-level series that extract_refit_summary.R omitted:
#   log10_conc_all_g   group-level wastewater concentration
#   mu_nb_g            group-level expected reported cases
#
# These are needed by the sliding-window correlation panel of manuscript
# figure 4 (panel B). Without them that panel is skipped and the panel letters
# shift.
#
# Run AFTER extract_refit_summary.R, on the same fit, from the project root.
#
# It uses the same thinning and the same chain ordering as
# extract_refit_summary.R, so the rows line up draw for draw with
# <prefix>_key_draws.rds and the two objects can be cbind-ed. If KEEP_EVERY
# below is changed, or the chains are combined in a different order, the rows
# will NOT correspond and the panel will plot mismatched draws without any
# error being raised.
###############################################################################

# ---- which fit to extract --------------------------------------------------
# Edit these two lines and run the script again for each fit.
#
#   "combined_fit_refit.RData"  prefix "refitted_combined"
#   "ww_fit_refit.RData"        prefix "refitted_ww"
#   "case_fit_refit.RData"      prefix "refitted_case"

fit_filename <- "combined_fit_refit.RData"
prefix       <- "refitted_combined"

# ----------------------------------------------------------------------------

fit_file <- file.path("outputs", "fits", fit_filename)
out_dir  <- file.path("outputs", "summaries")

KEEP_EVERY <- 20          # must match extract_refit_summary.R

extra <- c("log10_conc_all_g", "mu_nb_g")

if (!file.exists(fit_file)) stop("Cannot find ", fit_file, call. = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

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
message("  object: ", obj_name, "   chains: ", length(post),
        "   draws per chain: ", nrow(post[[1]]))

draws <- do.call(rbind, lapply(post, as.matrix))
message("  total draws: ", nrow(draws))

keep_cols <- unlist(lapply(extra, function(s)
  grep(paste0("^", s, "\\["), colnames(draws), value = TRUE)))

cat("\n--- nodes found ---\n")
for (s in extra) {
  n <- length(grep(paste0("^", s, "\\["), colnames(draws)))
  cat(sprintf("  %-18s %d columns\n", s, n))
  if (n == 0) {
    cat("    WARNING: not monitored in this fit. Figure 4 panel B cannot be\n")
    cat("    produced from it and the fit would need re-running with this node\n")
    cat("    added to the monitor list.\n")
  }
}

if (length(keep_cols) == 0) {
  stop("Neither series is present in this fit; nothing to write.", call. = FALSE)
}

rows  <- seq(1, nrow(draws), by = KEEP_EVERY)
small <- draws[rows, keep_cols, drop = FALSE]

# ---- row-alignment check against the key draws -----------------------------
# The figures script cbinds these two objects. If the row counts differ the
# cbind would recycle silently, so check here rather than discovering it in a
# plot.

key_file <- file.path(out_dir, paste0(prefix, "_key_draws.rds"))
if (file.exists(key_file)) {
  key <- readRDS(key_file)
  if (nrow(key) != nrow(small)) {
    stop("Row mismatch: ", basename(key_file), " has ", nrow(key),
         " draws but this object has ", nrow(small),
         ". KEEP_EVERY or the chain ordering must differ between the two ",
         "scripts. Do not use this file.", call. = FALSE)
  }
  message("  row alignment OK: ", nrow(small), " draws, matching ",
          basename(key_file))
} else {
  message("  note: ", basename(key_file),
          " not found, so row alignment was not checked")
}

saveRDS(small, file.path(out_dir, paste0(prefix, "_extra_draws.rds")))

cat("\n--- written ---\n")
cat("  ", file.path(out_dir, paste0(prefix, "_extra_draws.rds")),
    sprintf("  (%d draws x %d columns)", nrow(small), ncol(small)), "\n", sep = "")
cat(sprintf("     approximately %.1f MB\n",
            file.size(file.path(out_dir, paste0(prefix, "_extra_draws.rds"))) / 1e6))
