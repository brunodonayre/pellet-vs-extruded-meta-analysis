# Pelletized vs Extruded Feed Meta-analysis

This repository contains the R code used for the meta-analysis evaluating
the effect of feed processing technology (pelletization vs extrusion)
on specific growth rate (SGR) in aquaculture species.

## Software
- R (>= 4.2.0)
- Packages:
  - metafor
  - tidyverse
  - readxl

## Script
- `meta_analysis_sgr.R`: complete analysis pipeline, including:
  - data processing
  - effect size calculation (SMD)
  - random-effects meta-analysis
  - influence diagnostics
  - subgroup analyses (fish vs crustaceans)
  - publication bias assessment

## Data availability
The dataset used in this study is proprietary (Vita Science database)
and is not publicly available. Access can be granted upon reasonable request.

## Reproducibility
The script can be executed end-to-end once the dataset is placed
in the expected directory.

