###############################################################################
# derive_shedding_rates.R
#
# Purpose: derive the stage-specific viral shedding profile (phi_i) used in the
#          wastewater observation model, from measured viral loads reported in
#          Suner et al., Lancet Infect Dis 2023;23:445-453.
#
# Input:   suner_fig3a_digitised.csv, produced by digitise_suner_fig3a.py
# Output:  shedding_profile.csv  - the profile consumed by the fitting scripts
#
# ---------------------------------------------------------------------------
# WHY THIS SCRIPT EXISTS
#
# The previous derivation combined the two body sites by adding their log10
# concentrations:
#
#     phi_1 = 10^(7.3 + 5.0) = 10^12.3
#
# Adding logarithms multiplies the underlying concentrations rather than
# adding them. Two sites each shedding virus contribute an additive amount,
# so the concentrations must be summed on the linear scale:
#
#     phi_1 = 10^7.3 + 10^5.0 = 2.0e7   (log10 = 7.30)
#
# The previous derivation also took the stage-to-stage decline from Figure 2A
# (probability of detectable virus) rather than Figure 3A (measured viral
# load), and applied that probability to the log10 value rather than to the
# load. Figure 3A measures the quantity required, in bins that coincide with
# the model's three symptomatic stages, and is used here instead.
#
# ---------------------------------------------------------------------------
# UNITS: AN IMPORTANT LIMITATION
#
# Suner et al. report concentrations in log10 copies per mL OF VIRAL TRANSPORT
# MEDIA for lesion, pharynx and rectal swabs. These are not concentrations in
# stool or in lesion fluid, and they cannot be converted to copies per person
# per day without the swab recovery fraction, the transport medium volume, and
# the relationship between swab content and total daily output. None of these
# are available.
#
# The ABSOLUTE scale of phi_i is therefore not identifiable from this source.
# Only the RELATIVE profile across stages is.
#
# Accordingly this script returns phi normalised to stage 1, and the wastewater
# scaling parameter psi carries the entire absolute conversion. psi's prior in
# the fitting script must be widened correspondingly: the previous truncation
# spanned exactly one log10 unit, which was calibrated to the old phi scale and
# must not be carried over.
###############################################################################

# Base R only: no packages are required.

set.seed(20261102)

# ---- configuration ---------------------------------------------------------
# Paths follow the project layout: run this from the project root.

INPUT_CSV  <- file.path("outputs", "suner_fig3a_digitised.csv")
OUTPUT_CSV <- file.path("outputs", "shedding_profile.csv")

# Time bins of Figure 3A corresponding to the model's symptomatic stages.
# The model gives each of I1, I2, I3 a mean duration of roughly 7 days
# (theta_invall ~ dgamma(42.3, 1.99), mean 21.3 days, divided by three).
STAGE_BINS <- c(I1 = "1-7", I2 = "8-14", I3 = "15-21")

# Assay limit of detection for non-blood samples (Figure 3 caption).
LIMIT_OF_DETECTION <- 2.9

# Effective number of patients contributing to each binned median. Suner et al.
# followed 77 patients but do not report the number of samples per bin, so this
# is a conservative working value used only to set the width of the priors on
# the stage ratios. Vary it to test sensitivity.
N_EFFECTIVE <- 30

# Number of Monte Carlo draws used to propagate uncertainty into the ratios.
N_DRAWS <- 50000

# ---- helpers ---------------------------------------------------------------

#' Standard error of a median, from an interquartile range.
#'
#' For an approximately normal sample, SD ~ IQR / 1.349, and the standard error
#' of the median is about 1.253 * SD / sqrt(n).
median_se <- function(q1, q3, n) {
  sd_patient <- (q3 - q1) / 1.349
  1.253 * sd_patient / sqrt(n)
}

#' Combine site concentrations into one shedding quantity.
#'
#' Concentrations are added on the LINEAR scale. This is the step that was
#' previously performed on the log scale.
combine_sites <- function(log10_lesion, log10_rectum) {
  log10(10^log10_lesion + 10^log10_rectum)
}

# ---- load the digitised values ---------------------------------------------

if (!file.exists(INPUT_CSV)) {
  stop("Cannot find ", INPUT_CSV,
       ". Run digitise_suner_fig3a.py first.", call. = FALSE)
}

digitised <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)

#' Look up one site and bin, falling back to the limit of detection.
#'
#' Rectal medians fall to the limit of detection from bin 15-21 onwards, so no
#' box is drawn and nothing is extracted. Those observations are left censored:
#' the true median is at or below the limit of detection. Substituting the
#' limit itself is therefore an upper bound on that site's contribution.
lookup <- function(site, bin) {
  row <- digitised[digitised$site == site & digitised$bin == bin, ]
  if (nrow(row) == 0) {
    message("  note: ", site, " bin ", bin,
            " is at the limit of detection; using ", LIMIT_OF_DETECTION,
            " as an upper bound")
    return(list(median = LIMIT_OF_DETECTION,
                q1 = LIMIT_OF_DETECTION,
                q3 = LIMIT_OF_DETECTION,
                censored = TRUE))
  }
  list(median = row$median[1], q1 = row$q1[1], q3 = row$q3[1],
       censored = isTRUE(as.logical(row$censored[1])))
}

message("Reading digitised values from ", INPUT_CSV)
stages <- names(STAGE_BINS)
sites <- c("lesion", "rectum")

values <- list()
for (stage in stages) {
  for (site in sites) {
    values[[paste(site, stage, sep = "_")]] <- lookup(site, STAGE_BINS[[stage]])
  }
}

# ---- point estimates -------------------------------------------------------

phi_log10 <- vapply(stages, function(stage) {
  combine_sites(values[[paste0("lesion_", stage)]]$median,
                values[[paste0("rectum_", stage)]]$median)
}, numeric(1))

phi_relative <- 10^(phi_log10 - phi_log10[1])

# ---- uncertainty in the stage ratios ---------------------------------------
#
# The absolute level is carried by psi and is not estimated here, so only the
# ratios phi_2/phi_1 and phi_3/phi_1 need priors. Uncertainty in each binned
# median is propagated through the linear summation by Monte Carlo.

draws <- matrix(NA_real_, nrow = N_DRAWS, ncol = length(stages),
                dimnames = list(NULL, stages))

for (stage in stages) {
  lesion <- values[[paste0("lesion_", stage)]]
  rectum <- values[[paste0("rectum_", stage)]]

  lesion_draws <- rnorm(N_DRAWS, lesion$median,
                        median_se(lesion$q1, lesion$q3, N_EFFECTIVE))

  if (rectum$censored) {
    # Censored: the median lies at or below the limit of detection. Draw
    # uniformly over the interval below it, bounded well under the limit.
    rectum_draws <- runif(N_DRAWS, LIMIT_OF_DETECTION - 2, LIMIT_OF_DETECTION)
  } else {
    rectum_draws <- rnorm(N_DRAWS, rectum$median,
                          median_se(rectum$q1, rectum$q3, N_EFFECTIVE))
  }

  draws[, stage] <- combine_sites(lesion_draws, rectum_draws)
}

log10_ratio_draws <- draws[, stages[-1], drop = FALSE] - draws[, stages[1]]

ratio_prior <- data.frame(
  ratio = paste0("log10(phi_", 2:length(stages), "/phi_1)"),
  mean  = apply(log10_ratio_draws, 2, mean),
  sd    = apply(log10_ratio_draws, 2, sd),
  row.names = NULL
)

# ---- report ----------------------------------------------------------------

cat("\n")
cat("Stage-specific shedding, log10 copies per mL of viral transport media\n")
cat("---------------------------------------------------------------------\n")
cat(sprintf("%-6s %-8s %8s %8s %10s\n",
            "stage", "bin", "lesion", "rectum", "combined"))
for (i in seq_along(stages)) {
  stage <- stages[i]
  cat(sprintf("%-6s %-8s %8.2f %8.2f %10.3f\n",
              stage, STAGE_BINS[[stage]],
              values[[paste0("lesion_", stage)]]$median,
              values[[paste0("rectum_", stage)]]$median,
              phi_log10[i]))
}

cat("\nRelative profile, stage 1 = 1\n")
for (i in seq_along(stages)) {
  cat(sprintf("  phi_%d = %.6g\n", i, phi_relative[i]))
}

cat("\nStage ratios\n")
cat(sprintf("  phi_1/phi_2 = %.1fx\n", phi_relative[1] / phi_relative[2]))
cat(sprintf("  phi_2/phi_3 = %.0fx\n", phi_relative[2] / phi_relative[3]))
cat(sprintf("  phi_1/phi_3 = %.0fx  (total decline across the symptomatic period)\n",
            phi_relative[1] / phi_relative[3]))

cat("\nPriors for the stage ratios, if the profile is estimated rather than fixed\n")
cat(sprintf("  (N_EFFECTIVE = %d patients per binned median)\n", N_EFFECTIVE))
for (i in seq_len(nrow(ratio_prior))) {
  cat(sprintf("  %-22s  normal(mean = %7.4f, sd = %6.4f)\n",
              ratio_prior$ratio[i], ratio_prior$mean[i], ratio_prior$sd[i]))
}

cat("\nComparison with the values previously used\n")
previous <- c(12.3, 11.5, 10.2)
cat(sprintf("%-6s %12s %12s %10s\n", "stage", "previous", "corrected", "difference"))
for (i in seq_along(stages)) {
  cat(sprintf("%-6s %12s %12s %10s\n",
              stages[i],
              sprintf("10^%.2f", previous[i]),
              sprintf("10^%.2f", phi_log10[i]),
              sprintf("10^%.2f", previous[i] - phi_log10[i])))
}
cat(sprintf("\n  previous total decline across stages: %.0fx\n",
            10^(previous[1] - previous[3])))
cat(sprintf("  corrected total decline across stages: %.0fx\n",
            phi_relative[1] / phi_relative[3]))

# ---- write -----------------------------------------------------------------

output <- data.frame(
  stage          = stages,
  bin_days       = unname(STAGE_BINS[stages]),
  lesion_log10   = vapply(stages, function(s) values[[paste0("lesion_", s)]]$median, numeric(1)),
  rectum_log10   = vapply(stages, function(s) values[[paste0("rectum_", s)]]$median, numeric(1)),
  rectum_censored = vapply(stages, function(s) values[[paste0("rectum_", s)]]$censored, logical(1)),
  phi_log10      = unname(phi_log10),
  phi_relative   = unname(phi_relative),
  row.names = NULL
)

write.csv(output, OUTPUT_CSV, row.names = FALSE)
write.csv(ratio_prior, sub("\\.csv$", "_ratio_priors.csv", OUTPUT_CSV),
          row.names = FALSE)

cat("\nWritten:\n")
cat("  ", OUTPUT_CSV, "\n", sep = "")
cat("  ", sub("\\.csv$", "_ratio_priors.csv", OUTPUT_CSV), "\n", sep = "")

# ---- how to use this in the fitting script ---------------------------------

cat("\n")
cat("---------------------------------------------------------------------\n")
cat("To use this in fit_combined.R, replace the hardcoded block\n\n")
cat("    shed.I1 <- 10^12.3\n")
cat("    shed.I2 <- 10^11.5\n")
cat("    shed.I3 <- 10^10.2\n\n")
cat("with\n\n")
cat("    profile <- read.csv(file.path('outputs', 'shedding_profile.csv'))\n")
cat("    shed.I1 <- profile$phi_relative[1]\n")
cat("    shed.I2 <- profile$phi_relative[2]\n")
cat("    shed.I3 <- profile$phi_relative[3]\n\n")
cat("and widen the prior on the wastewater scaling parameter. The existing\n")
cat("prior\n\n")
cat("    log_mult ~ dnorm(log(3e-9), 2.5) T(log(1e-9), log(1e-8))\n\n")
cat("is truncated to exactly one log10 unit and was calibrated to the old phi\n")
cat("scale. With phi normalised to stage 1, psi must absorb the full\n")
cat("conversion from transport-medium concentration to copies per person per\n")
cat("day. A weakly informative replacement spanning several orders of\n")
cat("magnitude, for example\n\n")
cat("    log_mult ~ dnorm(log(3e3), 0.05) T(log(1), log(1e8))\n\n")
cat("leaves the scale to the data. Check the posterior does not sit against\n")
cat("either bound: that is what went undetected previously.\n")
cat("---------------------------------------------------------------------\n")
