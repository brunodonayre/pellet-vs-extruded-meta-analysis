# Pelletized vs. Extruded Feed Meta-analysis

This repository contains the **R code and derived summary outputs** used for the meta-analysis evaluating the effect of **feed processing technology (pelletization vs. extrusion)** on **specific growth rate (SGR)** in aquaculture species.

---

## Software

- **R ≥ 4.2.0**
- R packages:
  - `metafor`
  - `tidyverse`
  - `readxl`
  - `broom`
  - `ggplot2`

---

## Script

**`meta_analysis_sgr.R`**  
Complete, paper-ready analysis pipeline including:

- Data processing and filtering of accepted studies
- Effect size calculation with defined directionality  
  (Standardized Mean Difference, Hedges’ g, Log Response Ratio)
- Random-effects meta-analysis (REML)
- Sensitivity analyses (trimmed datasets)
- Multilevel models to account for non-independence
- Subgroup analyses (fish vs. crustaceans)
- Meta-regressions with biological and experimental moderators
- Publication bias assessment (funnel plots, Egger and Begg tests)
- Automated generation of supplementary tables and figures
- Risk-of-Bias assessment template

---

## Data availability

The dataset used in this study (**Vita Science database**) is **proprietary** and is **not publicly available**.

This repository does **not** contain raw experimental data.  
Only **derived summary tables and analytical outputs** are provided.

Access to the original dataset may be granted **upon reasonable request** to the corresponding authors.

---

## Reproducibility

The analysis script can be executed **end-to-end** once the proprietary dataset is placed in the expected directory.

All analytical decisions and model specifications are explicitly documented in the code, and package versions are recorded to ensure reproducibility.

