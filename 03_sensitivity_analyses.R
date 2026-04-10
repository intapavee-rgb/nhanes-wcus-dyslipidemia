# =============================================================================
# Script 03: Sensitivity Analyses (S1–S4)
# Project: Association between weekend catch-up sleep and dyslipidemia among U.S. adults: A cross-sectional study of NHANES 2017–March 2020
# Data: NHANES 2017–March 2020 Pre-Pandemic
# Author: Paveena Intaraksa
# =============================================================================

library(tidyverse)
library(survey)

# --- Load data and set up survey design ---
df <- read.csv("nhanes_analytic_dataset.csv")
df$WCUS_cat <- factor(df$WCUS_cat, levels = c("<0", "0to2", "2to4", "4plus"))

design <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA,
  weights = ~WTSAFPRP, nest = TRUE, data = df
)

# Full-model formula (used across all sensitivity analyses)
full_formula_cont <- dyslipidemia ~ WCUS + age + factor(sex) + factor(race) +
  factor(education) + PIR + factor(smoking) + factor(alcohol) +
  factor(phys_active) + weekday_sleep + BMI + factor(diabetes) + factor(hypertension)

full_formula_cat <- dyslipidemia ~ WCUS_cat + age + factor(sex) + factor(race) +
  factor(education) + PIR + factor(smoking) + factor(alcohol) +
  factor(phys_active) + weekday_sleep + BMI + factor(diabetes) + factor(hypertension)

# Helper: extract and print OR for continuous WCUS
print_cont_or <- function(model, label) {
  or <- round(exp(coef(model)["WCUS"]), 2)
  ci <- round(exp(confint(model)["WCUS", ]), 2)
  p  <- round(summary(model)$coefficients["WCUS", "Pr(>|t|)"], 3)
  cat(label, ": OR =", or, "(", ci[1], "–", ci[2], ") P =", p, "\n")
}

# Helper: extract and print ORs for categorical WCUS
print_cat_or <- function(model, label) {
  cat(label, ":\n")
  or <- round(exp(coef(model)), 2)
  ci <- round(exp(confint(model)), 2)
  for (cat_name in c("WCUS_cat0to2", "WCUS_cat2to4", "WCUS_cat4plus")) {
    p <- round(summary(model)$coefficients[cat_name, "Pr(>|t|)"], 3)
    cat("  ", cat_name, ": OR =", or[cat_name],
        "(", ci[cat_name, 1], "–", ci[cat_name, 2], ") P =", p, "\n")
  }
}

# =============================================================================
# S1: ALTERNATIVE REFERENCE CATEGORY (<=0 h only, excluding negative WCUS)
# =============================================================================
cat("=== S1: Alternative reference (<=0 h only) ===\n")

df_s1 <- df %>%
  mutate(WCUS_cat_s1 = factor(case_when(
    WCUS == 0              ~ "0",
    WCUS > 0 & WCUS < 2   ~ "0to2",
    WCUS >= 2 & WCUS < 4  ~ "2to4",
    WCUS >= 4              ~ "4plus",
    TRUE                   ~ NA_character_
  ), levels = c("0", "0to2", "2to4", "4plus"))) %>%
  filter(!is.na(WCUS_cat_s1))

design_s1 <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA,
  weights = ~WTSAFPRP, nest = TRUE, data = df_s1
)

m_s1 <- svyglm(
  dyslipidemia ~ WCUS_cat_s1 + age + factor(sex) + factor(race) +
    factor(education) + PIR + factor(smoking) + factor(alcohol) +
    factor(phys_active) + weekday_sleep + BMI + factor(diabetes) + factor(hypertension),
  design = design_s1, family = quasibinomial()
)
print_cat_or(m_s1, "S1 (ref = WCUS=0)")

# =============================================================================
# S2: EXCLUDING PARTICIPANTS ON LIPID-LOWERING MEDICATIONS
# =============================================================================
cat("\n=== S2: Excluding lipid-lowering medication users ===\n")

design_s2 <- subset(design, lipid_med == 0)
cat("S2 N =", nrow(design_s2$variables), "\n")

m_s2 <- svyglm(full_formula_cont, design = design_s2, family = quasibinomial())
print_cont_or(m_s2, "S2 (no lipid meds)")

# =============================================================================
# S3: EXCLUDING EXTREME WCUS VALUES (>= 6 h or <= -4 h)
# =============================================================================
cat("\n=== S3: Excluding extreme WCUS (>= 6 or <= -4) ===\n")

design_s3 <- subset(design, WCUS < 6 & WCUS > -4)
cat("S3 N =", nrow(design_s3$variables), "\n")

m_s3 <- svyglm(full_formula_cont, design = design_s3, family = quasibinomial())
print_cont_or(m_s3, "S3 (no extreme WCUS)")

# =============================================================================
# S4: MODEL WITHOUT DIABETES AND HYPERTENSION (POTENTIAL MEDIATORS)
# =============================================================================
cat("\n=== S4: Excluding diabetes and hypertension from model ===\n")

m_s4 <- svyglm(
  dyslipidemia ~ WCUS_cat + age + factor(sex) + factor(race) +
    factor(education) + PIR + factor(smoking) + factor(alcohol) +
    factor(phys_active) + weekday_sleep + BMI,
  design = design, family = quasibinomial()
)
print_cat_or(m_s4, "S4 (no diabetes/hypertension)")

cat("\nAll sensitivity analyses complete.\n")
sessionInfo()
