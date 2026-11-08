# ==================================================
# Script: derive_shedding_rates.R
# Purpose: Derive the stage-specific shedding rates (phi_i) used in the
#          wastewater observation model
# Project: Mpox_clean
# Updated: 2026-11-07
#
# Source: Suner et al. Viral dynamics in patients with monkeypox infection:
#         a prospective cohort study in Spain.
#         Lancet Infect Dis 2023;23:445-453, figure 3A.
# ==================================================

rm(list = ls())

options(scipen = 999)


# ---- Digitised values from figure 3A ----
# log10 copies per mL of viral transport medium, by body site and days since
# symptom onset. The three bins are the model's symptomatic stages I1, I2, I3.

fig3a <- read.csv(file.path("data", "suner_fig3a_digitised.csv"))

fig3a


# ---- Check against the values reported in the text of the paper ----
# Suner et al. report baseline medians of 7.3 (IQR 6.5-8.2) for lesion and
# 5.0 (IQR 2.9-7.5) for rectum.

check <- fig3a[fig3a$bin == "1-7", ]

cat("lesion baseline:", check$median[check$site == "lesion"],
    " paper reports 7.3\n")
cat("rectum baseline:", check$median[check$site == "rectum"],
    " paper reports 5.0\n")


# ---- Stage medians ----
lesion <- c(
  fig3a$median[fig3a$site == "lesion" & fig3a$bin == "1-7"],
  fig3a$median[fig3a$site == "lesion" & fig3a$bin == "8-14"],
  fig3a$median[fig3a$site == "lesion" & fig3a$bin == "15-21"]
)

# The rectal 15-21 bin has no box to read: its median sits at the assay limit
# of detection for non-blood samples (2.9, figure 3 caption). The true value
# is at or below this, so 2.9 is an upper bound on the rectal contribution.
limit_of_detection <- 2.9

rectum <- c(
  fig3a$median[fig3a$site == "rectum" & fig3a$bin == "1-7"],
  fig3a$median[fig3a$site == "rectum" & fig3a$bin == "8-14"],
  limit_of_detection
)


# ---- Combine the two sites on the linear scale ----
# The sites contribute additively to the virus shed by an infected individual,
# so the log10 values are converted to linear concentrations before summing.

phi <- 10^lesion + 10^rectum

shedding <- data.frame(
  stage       = c("I1", "I2", "I3"),
  bin         = c("1-7", "8-14", "15-21"),
  lesion      = round(lesion, 2),
  rectum      = round(rectum, 2),
  log10_phi   = round(log10(phi), 3),
  phi         = phi
)

shedding


# ---- Values used in the fitting scripts ----
cat("shed.I1 <- 10^", shedding$log10_phi[1], "\n", sep = "")
cat("shed.I2 <- 10^", shedding$log10_phi[2], "\n", sep = "")
cat("shed.I3 <- 10^", shedding$log10_phi[3], "\n", sep = "")

write.csv(shedding, file.path("outputs", "shedding_profile.csv"),
          row.names = FALSE)
