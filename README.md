# Pelletized vs. Extruded Feed Meta-analysis

This repository contains the **R code and derived summary outputs** used for the meta-analysis evaluating the effect of **feed processing technology (pelletization vs. extrusion)** on **specific growth rate (SGR)** in aquaculture species.

---

## Software

- **- R ≥ 4.2.0 (tested with R 4.5.1)**
- R packages:
  - `metafor`
  - `tidyverse`
  - `readxl`
  - `broom`
  - `ggplot2`

---

# Meta-analysis of SGR

## Repository structure

```text
.
├── meta_analysis_sgr.R
│   └── Main paper-ready R script implementing the full meta-analysis pipeline
│       (effect sizes, random-effects models, multilevel models, sensitivity
│       analyses, publication bias, and RoB template generation)
│
├── README.md
│   └── Project documentation and analysis overview
│
└── outputs/
    ├── Table_S2_sensitivity_models.csv
    │   └── Summary of sensitivity analyses (trimmed vs full datasets)
    │
    ├── Table_S3_Full_MetaRegression_Outputs.csv
    │   └── Complete meta-regression outputs for all tested moderators
    │
    ├── Table_S4_PublicationBias.csv
    │   └── Egger and Begg tests for publication bias
    │
    └── session_info.txt
        └── R session information (R version, OS, package versions)

```

## Note

Raw experimental data and comparison-level datasets are not included in this repository due to data-sharing restrictions.


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

## Analysis workflow (script outline)

The script `meta_analysis_sgr.R` is organized into the following sections:

1. **Setup**
   - Environment cleanup
   - Package loading
   - Definition of input/output paths

2. **Data processing**
   - Import of the proprietary Vita Science database
   - Filtering of accepted studies
   - Definition of effect size direction  
     *(positive values indicate an extrusion advantage)*

3. **Effect size calculation**
   - Standardized Mean Difference (SMD; Cohen’s d)
   - Hedges’ g (small-sample correction)
   - Log Response Ratio (LRR)

4. **Main meta-analyses**
   - Random-effects models (REML) for each effect size
   - Estimation of heterogeneity (τ², I²)
   - Prediction intervals

5. **Sensitivity analyses**
   - Exclusion of extreme comparisons (|SMD| > 5)
   - Re-estimation of main models on trimmed datasets

6. **Multilevel models**
   - Accounting for non-independence of multiple comparisons
   - Random effects nested by study and species
   - Subgroup analysis: fish vs. crustaceans

7. **Meta-regressions**
   - Biological and experimental moderators (e.g., species, duration, temperature, salinity)
   - Extraction of full regression outputs for supplementary tables

8. **Publication bias assessment**
   - Funnel plots
   - Egger’s regression test
   - Begg–Mazumdar rank correlation test

9. **Risk of Bias (RoB) framework**
   - Generation of a structured RoB assessment template
   - Optional visualization as a heatmap once ratings are completed

10. **Reproducibility**
    - Export of session information (`session_info.txt`)

---

## Outputs

The `outputs/` directory contains **derived analytical results only**, including:

- **Sensitivity model summaries** (Supplementary Table S2)
- **Full meta-regression outputs** (Supplementary Table S3)
- **Publication bias statistics** (Supplementary Table S4)
- **Session information** documenting the computational environment

These files correspond to supplementary materials used in the associated manuscript.

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

