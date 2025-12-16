############################################################
## Pelletized vs Extruded Feed – Meta-analysis on SGR
## Clean / Paper-ready analysis pipeline
## - Direction: positive effect = extrusion advantage
## - Effect sizes: SMD (Cohen's d), Hedges' g, Log Response Ratio (LRR)
## - Models: Random-effects (REML), Multilevel (rma.mv)
## - Sensitivity: trimmed |SMD| < 5
## - Outputs: CSV tables + Funnel plots + (optional) RoB template
############################################################

##### ======================================================
##### 1) Setup
##### ======================================================

rm(list = ls())
graphics.off()

library(readxl)
library(tidyverse)
library(metafor)
library(broom)
library(ggplot2)

# ---- Paths ----
PATH_DB  <- file.path("data", "Vita_science.xlsx")
SHEET_DB <- "db_final"


# Output folder (will be created if not exists)
OUT_DIR <- "outputs_meta_analysis"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)


if (!file.exists(PATH_DB)) {
  stop(
    "Data file not found: ", PATH_DB,
    "\nPlease place 'Vita_science.xlsx' inside the 'data/' folder."
  )
}


##### ======================================================
##### 2) Data processing
##### ======================================================



db_raw <- read_excel(path = PATH_DB, sheet = SHEET_DB)

db <- db_raw %>%
  filter(rechazado_acceptado == "aceptado")

# Effect size direction:
#   Positive effect => Extrusion advantage
#   Group1 = Extruded, Group2 = Pelletized
db <- escalc(
  measure = "SMD",
  n1i  = n_ind_ext,  n2i  = n_ind_pel,
  m1i  = mean_sgr_ext, m2i = mean_sgr_pel,
  sd1i = std_sgr_ext,  sd2i = std_sgr_pel,
  data = db
)

# Add sensitivity metrics:
# - Hedges' g (small-sample correction from yi/vi)
# - LRR = log(mean_ext / mean_pel) + variance
db <- db %>%
  mutate(
    J = 1 - (3 / (4 * (n_ind_ext + n_ind_pel - 2) - 1)),
    hedges_g = yi * J,
    var_g    = vi * J^2,
    lrr      = log(mean_sgr_ext / mean_sgr_pel),
    var_lrr  = (std_sgr_ext^2) / (n_ind_ext * mean_sgr_ext^2) +
      (std_sgr_pel^2) / (n_ind_pel * mean_sgr_pel^2)
  )

# Quick check
summary(as.data.frame(db)[, c("yi", "hedges_g", "lrr")])

##### ======================================================
##### 3) Main random-effects meta-analyses (Full dataset)
##### ======================================================

res_smd <- rma(yi = yi,       vi = vi,      data = db, method = "REML")
res_g   <- rma(yi = hedges_g, vi = var_g,   data = db, method = "REML")
res_lrr <- rma(yi = lrr,      vi = var_lrr, data = db, method = "REML")

cat("\n--- FULL DATASET MODELS ---\n")
print(summary(res_smd))
print(summary(res_g))
print(summary(res_lrr))

# Prediction intervals (commonly reported)
cat("\n--- PREDICTION INTERVALS (Hedges' g & LRR) ---\n")
print(predict(res_g))
print(predict(res_lrr))

##### ======================================================
##### 4) Export: included comparisons (Supplementary Table S1)
##### ======================================================

table_S1 <- db %>%
  select(study_id, doi, especie, especie_grupo,
         mean_sgr_ext, std_sgr_ext, n_ind_ext,
         mean_sgr_pel, std_sgr_pel, n_ind_pel,
         yi, vi, hedges_g, var_g, lrr, var_lrr)

write.csv(table_S1, file.path(OUT_DIR, "Table_S1_included_comparisons.csv"), row.names = FALSE)

# Summary of models (Supplementary Table S2)
model_summary <- data.frame(
  Metric   = c("SMD", "Hedges_g", "LRR"),
  Estimate = c(as.numeric(res_smd$b), as.numeric(res_g$b), as.numeric(res_lrr$b)),
  CI_Lower = c(res_smd$ci.lb, res_g$ci.lb, res_lrr$ci.lb),
  CI_Upper = c(res_smd$ci.ub, res_g$ci.ub, res_lrr$ci.ub),
  Tau2     = c(res_smd$tau2,  res_g$tau2,  res_lrr$tau2),
  I2       = c(res_smd$I2,    res_g$I2,    res_lrr$I2)
)

write.csv(model_summary, file.path(OUT_DIR, "Table_S2_sensitivity_models.csv"), row.names = FALSE)

##### ======================================================
##### 5) Sensitivity analysis: trim extreme SMDs (|SMD| < 5)
##### ======================================================

# Identify extreme comparisons (informative)
extreme_rows <- db %>%
  filter(yi > 5 | yi < -5) %>%
  select(study_id, yi, mean_sgr_ext, mean_sgr_pel)

if (nrow(extreme_rows) > 0) {
  write.csv(extreme_rows, file.path(OUT_DIR, "extreme_comparisons_absSMD_gt5.csv"), row.names = FALSE)
}

db_trim <- db %>% filter(abs(yi) < 5)

res_g_trim   <- rma(yi = hedges_g, vi = var_g,   data = db_trim, method = "REML")
res_lrr_trim <- rma(yi = lrr,      vi = var_lrr, data = db_trim, method = "REML")

cat("\n--- TRIMMED DATASET MODELS (|SMD|<5) ---\n")
print(summary(res_g_trim))
print(predict(res_g_trim))
print(summary(res_lrr_trim))
print(predict(res_lrr_trim))

##### ======================================================
##### 6) Multilevel model (non-independence) on trimmed dataset
##### ======================================================

# Multilevel: random effects within study and species (nested)
# Moderator: species group (Fish vs Crustacean)
res_ml <- rma.mv(
  yi = hedges_g,
  V  = var_g,
  random = ~ 1 | study_id/especie,
  mods   = ~ especie_grupo,
  data   = db_trim,
  method = "REML"
)

cat("\n--- MULTILEVEL MODEL (trimmed) ---\n")
print(summary(res_ml))

# Predicted marginal means for groups (safe construction of newmods)
# For mods = ~ especie_grupo, newmods should have the columns of the model matrix without intercept.
mm <- model.matrix(~ especie_grupo, data = data.frame(especie_grupo = c("Crustacean", "Fish")))
newmods <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]

pred_groups <- predict(res_ml, newmods = newmods)

pred_df <- data.frame(
  especie_grupo = c("Crustacean", "Fish"),
  estimate = pred_groups$pred,
  ci.lb    = pred_groups$ci.lb,
  ci.ub    = pred_groups$ci.ub
)

# Plot: jitter + group means (multilevel)
p_group <- ggplot() +
  geom_jitter(
    data = db_trim,
    aes(x = especie_grupo, y = hedges_g),
    width = 0.15, alpha = 0.5, size = 2
  ) +
  geom_pointrange(
    data = pred_df,
    aes(x = especie_grupo, y = estimate, ymin = ci.lb, ymax = ci.ub),
    size = 1.1, fatten = 2.5
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Species group",
    y = "Effect size (Hedges' g)",
    title = "Extruded vs. pelleted feeds: Fish vs. Crustaceans",
    subtitle = "Points: individual comparisons; Lines: multilevel model means ±95% CI"
  ) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "Figure_species_group_multilevel.png"), p_group, width = 10, height = 7, dpi = 300)

##### ======================================================
##### 7) Moderators (meta-regressions) on trimmed dataset
##### ======================================================

# NOTE: These require that the columns exist and have usable variation.
# If a moderator is missing, comment that model out.

res_year <- rma(yi = hedges_g, vi = var_g, mods = ~ ano, data = db_trim, method = "REML")
res_temp <- rma(yi = hedges_g, vi = var_g, mods = ~ as.factor(temperatura), data = db_trim, method = "REML")
res_sal  <- rma(yi = hedges_g, vi = var_g, mods = ~ as.factor(salinidad),   data = db_trim, method = "REML")
res_dur  <- rma(yi = hedges_g, vi = var_g, mods = ~ duracion,               data = db_trim, method = "REML")
res_species <- rma(yi = hedges_g, vi = var_g, mods = ~ especie, data = db_trim, method = "REML")

# Optional: random by study/species (no fixed moderator here, just variance structure)
res_ml_species <- rma.mv(
  yi = hedges_g, V = var_g,
  random = ~ 1 | study_id/as.factor(especie),
  data = db_trim, method = "REML"
)

# Helper to extract results consistently
extract_mr <- function(model, model_name){
  
  td <- broom::tidy(model, conf.int = TRUE)
  
  if(!("conf.low" %in% names(td))){
    ci <- confint(model)
    td$conf.low  <- ci$ci.lb
    td$conf.high <- ci$ci.ub
  }
  
  tau2_val <- NA
  if(!is.null(model$tau2))   tau2_val <- model$tau2
  if(!is.null(model$sigma2)) tau2_val <- model$sigma2[1]
  
  I2_val <- ifelse(!is.null(model$I2), model$I2, NA)
  R2_val <- ifelse(!is.null(model$R2), model$R2, NA)
  
  td$Model <- model_name
  td$tau2  <- tau2_val
  td$I2    <- I2_val
  td$R2    <- R2_val
  
  td %>%
    dplyr::select(Model, term, estimate, std.error, statistic, p.value,
                  conf.low, conf.high, tau2, I2, R2)
}

tab_year    <- extract_mr(res_year,    "Year")
tab_temp    <- extract_mr(res_temp,    "Temperature")
tab_sal     <- extract_mr(res_sal,     "Salinity")
tab_dur     <- extract_mr(res_dur,     "Duration")
tab_species <- extract_mr(res_species, "Species")
tab_ml_grp  <- extract_mr(res_ml,      "Species_group_ML")
tab_ml_sp   <- extract_mr(res_ml_species, "Random_study_species_ML")

Table_S3 <- bind_rows(tab_year, tab_temp, tab_sal, tab_dur, tab_species, tab_ml_grp, tab_ml_sp) %>%
  rename(
    Term = term,
    Estimate = estimate,
    SE = std.error,
    Z = statistic,
    P = p.value,
    CI_Lower = conf.low,
    CI_Upper = conf.high
  )

write.csv(Table_S3, file.path(OUT_DIR, "Table_S3_Full_MetaRegression_Outputs.csv"), row.names = FALSE)

##### ======================================================
##### 8) Publication bias (Full vs Trimmed) using Hedges' g models
##### ======================================================

# Funnel plots
png(file.path(OUT_DIR, "Funnel_FullDataset.png"), width = 900, height = 900, res = 150)
funnel(res_g, main = paste0("Funnel plot – Full dataset (Hedges' g), k = ", res_g$k))
dev.off()

png(file.path(OUT_DIR, "Funnel_TrimmedDataset.png"), width = 900, height = 900, res = 150)
funnel(res_g_trim, main = paste0("Funnel plot – Trimmed dataset (Hedges' g), k = ", res_g_trim$k))
dev.off()

# Egger + Begg
egger_full <- regtest(res_g,     model = "rma", predictor = "sei")
begg_full  <- ranktest(res_g)

egger_trim <- regtest(res_g_trim, model = "rma", predictor = "sei")
begg_trim  <- ranktest(res_g_trim)

Table_S4 <- data.frame(
  Dataset  = c("Full", "Full", "Trimmed", "Trimmed"),
  Test     = c("Egger", "Begg", "Egger", "Begg"),
  Statistic = c(egger_full$zval, begg_full$tau, egger_trim$zval, begg_trim$tau),
  P_value   = c(egger_full$pval, begg_full$pval, egger_trim$pval, begg_trim$pval)
)

write.csv(Table_S4, file.path(OUT_DIR, "Table_S4_PublicationBias.csv"), row.names = FALSE)

cat("\n--- PUBLICATION BIAS SUMMARY ---\n")
cat("Full dataset (k =", res_g$k, "): Egger z =", round(egger_full$zval,3), "p =", signif(egger_full$pval,3),
    " | Begg tau =", round(begg_full$tau,3), "p =", signif(begg_full$pval,3), "\n")
cat("Trimmed dataset (k =", res_g_trim$k, "): Egger z =", round(egger_trim$zval,3), "p =", signif(egger_trim$pval,3),
    " | Begg tau =", round(begg_trim$tau,3), "p =", signif(begg_trim$pval,3), "\n")

##### ======================================================
##### 9) Risk of Bias (RoB) template export (manual filling)
##### ======================================================

study_list <- unique(db$study_id)

rob_domains <- c(
  "Randomization_of_tanks",
  "Allocation_concealment",
  "Blinding_of_feeders/assessors",
  "Baseline_equivalence",
  "Incomplete_outcome_data",
  "Selective_reporting",
  "Feed_preparation_reporting",
  "Water_quality_control",
  "Tank_replication_and_pseudoreplication"
)

Table_S5_RoB <- expand.grid(
  Study  = study_list,
  Domain = rob_domains
)

# Fill manually later with: Low / Unclear / High
Table_S5_RoB$Rating <- ""

write.csv(Table_S5_RoB, file.path(OUT_DIR, "Table_S5_RiskOfBias_Template.csv"), row.names = FALSE)

# Optional: only plot heatmap if you have filled ratings
if (any(Table_S5_RoB$Rating %in% c("Low","Unclear","High"))) {
  
  rob_plot <- Table_S5_RoB %>%
    mutate(
      Rating = factor(Rating, levels = c("Low", "Unclear", "High"))
    ) %>%
    ggplot(aes(x = Domain, y = Study, fill = Rating)) +
    geom_tile(color = "grey70") +
    scale_fill_manual(values = c("Low" = "forestgreen",
                                 "Unclear" = "gold",
                                 "High" = "firebrick"),
                      drop = FALSE) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(title = "Risk of Bias Across Included Studies",
         x = "Risk-of-Bias Domain",
         y = "Study ID")
  
  ggsave(file.path(OUT_DIR, "Figure_RoB_Heatmap.png"), rob_plot, width = 12, height = 8, dpi = 300)
  
} else {
  message("RoB heatmap not generated because 'Rating' is empty. Fill Table_S5_RiskOfBias_Template.csv and re-run the RoB plotting block.")
}

##### ======================================================
##### 10) Reproducibility (optional)
##### ======================================================

sink(file.path(OUT_DIR, "session_info.txt"))
sessionInfo()
sink()

cat("\nDONE. Outputs saved in: ", OUT_DIR, "\n")
