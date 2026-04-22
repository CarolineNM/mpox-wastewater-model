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

- `scripts/` – Model fitting, diagnostics, and analysis scripts  
- `data/` – Placeholder directory (raw surveillance data not included)  
- `outputs/` – Model outputs (not version controlled)  

## Reproducibility Workflow

To reproduce the full analysis:

1. `scripts/fit_case_only.R`
2. `scripts/fit_ww_only.R`
3. `scripts/fit_combined.R`
4. `scripts/generate_manuscript_figures.R`
5. `scripts/run_counterfactuals.R`

## Data Availability

Raw surveillance and wastewater data are not included due to data governance restrictions.  
Processed data may be provided upon reasonable request.

## Software Requirements

- R (≥ 4.4.0)
- JAGS
- Required R packages listed within scripts