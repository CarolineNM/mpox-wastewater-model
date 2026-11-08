# Mpox Wastewater Transmission Model

This repository contains code supporting the manuscript:

**Bridging wastewater and case surveillance reveals the hidden
burden of the 2022 mpox outbreak**

## Overview

We developed a deterministic compartmental transmission model integrating:

- Reported case surveillance data  
- Wastewater viral load measurements  
- Vaccination dynamics  

The model was fitted using Bayesian inference in JAGS and implemented in R.

## Repository Structure

- `scripts/` – Model fitting, extraction and analysis scripts  
- `data/` – Digitised shedding values; raw surveillance data not included  
- `outputs/` – Model outputs (not version controlled)  

## Reproducibility Workflow

Run from the project root, in this order.

**1. Shedding rates**

- `derive_shedding_rates.R` – derives the stage-specific shedding rates from
  the digitised values in `data/suner_fig3a_digitised.csv`

**2. Model fitting**

- `fit_case_only.R`, `fit_ww_only.R`, `fit_combined.R` – the three model variants
- `fit_combined_low.R`, `fit_combined_high.R` – viral decay sensitivity
- `fit_combined_loshed.R`, `fit_combined_highshed.R`, `fit_combined_flatshed.R` –
  vaccinated shedding sensitivity

**3. Extraction**

- `extract_refit_summary.R` – posterior summaries, convergence diagnostics
  and thinned draws, run once per fit
- `extract_extra_series.R` – group-level series for figure 4B
- `extract_sensitivity.R` – burden estimates across the sensitivity fits,
  for tables S2 and S3

**4. Analysis and figures**

- `analyse_refit.R` – burden estimates and underascertainment factors
- `generate_manuscript_figures.R` – manuscript and supplementary figures
- `run_counterfactuals.R` – counterfactual scenario analyses
- `sero_validation.R` – comparison with observed seroprevalence
- `Sens_run_counterfactuals.R`, `decay_rateimpact.R`, `shedding_impact.R` –
  sensitivity analyses

## Data availability

The mpox case-based surveillance, wastewater and vaccination datasets are not
openly available and may be available from the corresponding author on
reasonable request. They are held in controlled-access storage at the British
Columbia Centre for Disease Control.

`data/suner_fig3a_digitised.csv` holds the viral load values extracted from
figure 3A of Suner et al. (Lancet Infect Dis 2023;23:445-453), used to derive
the stage-specific shedding rates.
