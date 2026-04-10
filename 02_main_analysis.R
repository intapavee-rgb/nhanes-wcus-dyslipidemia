# =============================================================================
# Script 02: Main Analysis
# Project: Association between weekend catch-up sleep and dyslipidemia among U.S. adults: A cross-sectional study of NHANES 2017–March 2020
# Data: NHANES 2017–March 2020 Pre-Pandemic
# Author: Paveena Intaraksa
# =============================================================================

library(tidyverse)
library(survey)

# =============================================================================
# SECTION 1: LOAD DATA AND SET UP SURVEY DESIGN
# =============================================================================

df <- read.csv("nhanes_analytic_dataset.csv")

# Set factor levels
df$WCUS_cat <- factor(df$WCUS_cat, levels = c("<0", "0to2", "2to4", "4plus"))

# Complex survey design using fasting subsample weights
design <- svydesign(
  id      = ~SDMVPSU,    # Primary sampling unit
  strata  = ~SDMVSTRA,   # Stratification variable
  weights = ~WTSAFPRP,   # Fasting subsample weight
  nest    = TRUE,
  data    = df
)

cat("Analytic N =", nrow(df), "\n")

# =============================================================================
# SECTION 2: TABLE 1 — WEIGHTED BASELINE CHARACTERISTICS
# =============================================================================

cat("\n--- TABLE 1: BASELINE CHARACTERISTICS ---\n")

# Overall weighted N
cat("Weighted N:", round(sum(weights(design)), 0), "\n")

# Weighted dyslipidemia prevalence
dys_prev <- svymean(~dyslipidemia, design)
cat("Dyslipidemia prevalence:", round(coef(dys_prev) * 100, 1), "% (SE",
    round(SE(dys_prev) * 100, 1), "%)\n")

# Age
age_est <- svymean(~age, design)
cat("Mean age:", round(coef(age_est), 1), "(SE", round(SE(age_est), 1), ")\n")

# Sex (1 = Male)
sex_est <- svymean(~I(sex == 1), design)
cat("Male %:", round(coef(sex_est) * 100, 1), "\n")

# WCUS categories
wcus_table <- svytable(~WCUS_cat, design)
cat("\nWCUS distribution:\n")
print(round(prop.table(wcus_table) * 100, 1))

# =============================================================================
# SECTION 3: PRIMARY LOGISTIC REGRESSION — MODELS 1–3 (TABLE 2)
# =============================================================================

cat("\n--- TABLE 2: LOGISTIC REGRESSION MODELS ---\n")

# Model 1: Demographic adjustment
m1 <- svyglm(
  dyslipidemia ~ WCUS_cat + age + factor(sex) + factor(race) +
    factor(education) + PIR,
  design = design,
  family = quasibinomial()
)

# Model 2: + Lifestyle factors
m2 <- svyglm(
  dyslipidemia ~ WCUS_cat + age + factor(sex) + factor(race) +
    factor(education) + PIR +
    factor(smoking) + factor(alcohol) + factor(phys_active) + weekday_sleep,
  design = design,
  family = quasibinomial()
)

# Model 3: Full adjustment (primary model)
m3 <- svyglm(
  dyslipidemia ~ WCUS_cat + age + factor(sex) + factor(race) +
    factor(education) + PIR +
    factor(smoking) + factor(alcohol) + factor(phys_active) + weekday_sleep +
    BMI + factor(diabetes) + factor(hypertension),
  design = design,
  family = quasibinomial()
)

# Extract and display ORs with 95% CIs
extract_or <- function(model, label) {
  cat("\n", label, ":\n", sep = "")
  or  <- round(exp(coef(model)), 2)
  ci  <- round(exp(confint(model)), 2)
  for (cat_name in c("WCUS_cat0to2", "WCUS_cat2to4", "WCUS_cat4plus")) {
    p <- round(summary(model)$coefficients[cat_name, "Pr(>|t|)"], 3)
    cat("  ", cat_name, ": OR =", or[cat_name],
        "(", ci[cat_name, 1], "–", ci[cat_name, 2], ") P =", p, "\n")
  }
}

extract_or(m1, "Model 1")
extract_or(m2, "Model 2")
extract_or(m3, "Model 3")

# --- Continuous WCUS (Model 3) ---
m3_cont <- svyglm(
  dyslipidemia ~ WCUS + age + factor(sex) + factor(race) +
    factor(education) + PIR +
    factor(smoking) + factor(alcohol) + factor(phys_active) + weekday_sleep +
    BMI + factor(diabetes) + factor(hypertension),
  design = design,
  family = quasibinomial()
)
or_cont <- round(exp(coef(m3_cont)["WCUS"]), 2)
ci_cont <- round(exp(confint(m3_cont)["WCUS", ]), 2)
p_cont  <- round(summary(m3_cont)$coefficients["WCUS", "Pr(>|t|)"], 3)
cat("\nContinuous WCUS (Model 3): OR =", or_cont,
    "(", ci_cont[1], "–", ci_cont[2], ") P =", p_cont, "\n")

# --- P for trend (ordinal WCUS) ---
m3_trend <- svyglm(
  dyslipidemia ~ WCUS_ordinal + age + factor(sex) + factor(race) +
    factor(education) + PIR +
    factor(smoking) + factor(alcohol) + factor(phys_active) + weekday_sleep +
    BMI + factor(diabetes) + factor(hypertension),
  design = design,
  family = quasibinomial()
)
p_trend <- round(summary(m3_trend)$coefficients["WCUS_ordinal", "Pr(>|t|)"], 3)
cat("P for trend:", p_trend, "\n")

# =============================================================================
# SECTION 4: STRATIFIED ANALYSES — TABLE 3 (CONTINUOUS WCUS)
# =============================================================================

cat("\n--- TABLE 3: STRATIFIED ANALYSES (CONTINUOUS WCUS) ---\n")

run_stratified_cont <- function(subset_condition, label) {
  d_sub <- subset(design, eval(parse(text = subset_condition)))
  m <- svyglm(
    dyslipidemia ~ WCUS + age + factor(sex) + factor(race) +
      factor(education) + PIR +
      factor(smoking) + factor(alcohol) + factor(phys_active) + weekday_sleep +
      BMI + factor(diabetes) + factor(hypertension),
    design = d_sub, family = quasibinomial()
  )
  or <- round(exp(coef(m)["WCUS"]), 2)
  ci <- round(exp(confint(m)["WCUS", ]), 2)
  p  <- round(summary(m)$coefficients["WCUS", "Pr(>|t|)"], 3)
  cat(label, ": OR =", or, "(", ci[1], "–", ci[2], ") P =", p, "\n")
}

# By sex
run_stratified_cont("sex == 1", "Male")
run_stratified_cont("sex == 2", "Female")

# By weekday sleep (short = <7 h, adequate = >=7 h)
run_stratified_cont("weekday_sleep < 7",  "Short weekday sleep (<7 h)")
run_stratified_cont("weekday_sleep >= 7", "Adequate weekday sleep (≥7 h)")

# By age group
run_stratified_cont("age < 45",  "Age <45")
run_stratified_cont("age >= 45", "Age ≥45")

# By BMI group
run_stratified_cont("BMI < 25",  "Normal weight (BMI <25)")
run_stratified_cont("BMI >= 25", "Overweight/obese (BMI ≥25)")

# =============================================================================
# SECTION 5: SEX-STRATIFIED CATEGORICAL — TABLE 4
# =============================================================================

cat("\n--- TABLE 4: SEX-STRATIFIED CATEGORICAL WCUS ---\n")

run_stratified_cat <- function(sex_val, label) {
  d_sub <- subset(design, sex == sex_val)
  m <- svyglm(
    dyslipidemia ~ WCUS_cat + age + factor(race) +
      factor(education) + PIR +
      factor(smoking) + factor(alcohol) + factor(phys_active) + weekday_sleep +
      BMI + factor(diabetes) + factor(hypertension),
    design = d_sub, family = quasibinomial()
  )
  cat("\n", label, ":\n")
  or <- round(exp(coef(m)), 2)
  ci <- round(exp(confint(m)), 2)
  for (cat_name in c("WCUS_cat0to2", "WCUS_cat2to4", "WCUS_cat4plus")) {
    p <- round(summary(m)$coefficients[cat_name, "Pr(>|t|)"], 3)
    cat("  ", cat_name, ": OR =", or[cat_name],
        "(", ci[cat_name, 1], "–", ci[cat_name, 2], ") P =", p, "\n")
  }
}

run_stratified_cat(1, "Males")
run_stratified_cat(2, "Females")

# =============================================================================
# SECTION 6: INDIVIDUAL LIPID COMPONENTS — TABLE 5
# =============================================================================

cat("\n--- TABLE 5: INDIVIDUAL LIPID COMPONENTS ---\n")

# Linear regression for each continuous lipid outcome
lipid_outcomes <- c("LBXTC", "LBDLDL", "LBDHDD", "LBXTR")
lipid_labels   <- c("Total cholesterol", "LDL-C", "HDL-C", "Triglycerides")

for (i in seq_along(lipid_outcomes)) {
  outcome <- lipid_outcomes[i]
  if (!outcome %in% names(df)) next
  m <- svyglm(
    as.formula(paste0(
      outcome, " ~ WCUS + age + factor(sex) + factor(race) + factor(education) + PIR +",
      "factor(smoking) + factor(alcohol) + factor(phys_active) + weekday_sleep +",
      "BMI + factor(diabetes) + factor(hypertension)"
    )),
    design = design,
    family = gaussian()
  )
  beta <- round(coef(m)["WCUS"], 3)
  ci   <- round(confint(m)["WCUS", ], 3)
  p    <- round(summary(m)$coefficients["WCUS", "Pr(>|t|)"], 3)
  cat(lipid_labels[i], ": Beta =", beta,
      "(", ci[1], "–", ci[2], ") P =", p, "\n")
}

# Session info
sessionInfo()
