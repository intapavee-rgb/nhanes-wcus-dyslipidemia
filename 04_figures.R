# =============================================================================
# Script 04: Figure Generation
# Project: Association between weekend catch-up sleep and dyslipidemia among U.S. adults: A cross-sectional study of NHANES 2017–March 2020
# Data: NHANES 2017–March 2020 Pre-Pandemic
# Author: Paveena Intaraksa
# =============================================================================

# install.packages(c("survey", "rms", "ggplot2", "tidyverse"))

library(tidyverse)
library(survey)
library(rms)
library(ggplot2)

# --- Load data and survey design ---
df <- read.csv("nhanes_analytic_dataset.csv")
df$WCUS_cat <- factor(df$WCUS_cat, levels = c("<0", "0to2", "2to4", "4plus"))

design <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA,
  weights = ~WTSAFPRP, nest = TRUE, data = df
)

# =============================================================================
# FIGURE 1: PARTICIPANT FLOW DIAGRAM
# =============================================================================
# Figure 1 was created manually as a PowerPoint/Word diagram.
# Key numbers:
#   NHANES 2017-March 2020 examined: [total]
#   Excluded age <18: [n]
#   Excluded pregnant: [n]
#   Excluded missing WCUS: [n]
#   Excluded missing lipid data: [n]
#   Excluded missing covariates (primarily PIR): ~805
#   Final analytic sample: 2,979

# =============================================================================
# FIGURE 2: RESTRICTED CUBIC SPLINE (RCS) — DOSE-RESPONSE
# =============================================================================

cat("Generating Figure 2: RCS dose-response plot...\n")

# Set up datadist for rms package
dd <- datadist(df)
options(datadist = "dd")

# Fit RCS model with 4 knots using survey-weighted approach
# Note: rms does not directly support svydesign; use weighted GLM as approximation
rcs_model <- glm(
  dyslipidemia ~ rcs(WCUS, 4) + age + factor(sex) + factor(race) +
    factor(education) + PIR + factor(smoking) + factor(alcohol) +
    factor(phys_active) + weekday_sleep + BMI + factor(diabetes) + factor(hypertension),
  data    = df,
  weights = WTSAFPRP,
  family  = quasibinomial()
)

# Prediction grid
wcus_grid <- seq(min(df$WCUS, na.rm = TRUE), max(df$WCUS, na.rm = TRUE), length.out = 200)

new_data <- data.frame(
  WCUS          = wcus_grid,
  age           = mean(df$age, na.rm = TRUE),
  sex           = 1,
  race          = 1,
  education     = 3,
  PIR           = mean(df$PIR, na.rm = TRUE),
  smoking       = 0,
  alcohol       = 1,
  phys_active   = 1,
  weekday_sleep = mean(df$weekday_sleep, na.rm = TRUE),
  BMI           = mean(df$BMI, na.rm = TRUE),
  diabetes      = 0,
  hypertension  = 0
)

pred <- predict(rcs_model, newdata = new_data, type = "link", se.fit = TRUE)
new_data$log_or  <- pred$fit
new_data$se      <- pred$se.fit
new_data$OR      <- exp(new_data$log_or)
new_data$OR_low  <- exp(new_data$log_or - 1.96 * new_data$se)
new_data$OR_high <- exp(new_data$log_or + 1.96 * new_data$se)

# Reference: OR = 1 at WCUS = 0
ref_log_or <- new_data$log_or[which.min(abs(new_data$WCUS))]
new_data$OR      <- exp(new_data$log_or - ref_log_or)
new_data$OR_low  <- exp(new_data$log_or - ref_log_or - 1.96 * new_data$se)
new_data$OR_high <- exp(new_data$log_or - ref_log_or + 1.96 * new_data$se)

# WCUS histogram for rug/density
wcus_hist <- df %>%
  count(WCUS = round(WCUS)) %>%
  mutate(pct = n / sum(n) * 100)

# Plot
fig2 <- ggplot(new_data, aes(x = WCUS, y = OR)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = OR_low, ymax = OR_high), fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue", linewidth = 1) +
  scale_x_continuous(breaks = seq(-4, 8, by = 2),
                     limits = c(-4, 8)) +
  scale_y_continuous(limits = c(0.5, 2.5)) +
  labs(
    x = "Weekend Catch-Up Sleep Duration (hours)",
    y = "Odds Ratio for Dyslipidemia (95% CI)",
    title = "Dose-Response Relationship between WCUS and Dyslipidemia"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 11)
  )

ggsave("Figure2_RCS.pdf", fig2, width = 7, height = 5, dpi = 300)
ggsave("Figure2_RCS.png", fig2, width = 7, height = 5, dpi = 300)
cat("Saved: Figure2_RCS.pdf / .png\n")

# =============================================================================
# FIGURE 3: FOREST PLOT — WCUS CATEGORIES ACROSS MODELS 1–3
# =============================================================================

cat("Generating Figure 3: Forest plot...\n")

# Fit models
m1 <- svyglm(
  dyslipidemia ~ WCUS_cat + age + factor(sex) + factor(race) + factor(education) + PIR,
  design = design, family = quasibinomial()
)
m2 <- svyglm(
  dyslipidemia ~ WCUS_cat + age + factor(sex) + factor(race) + factor(education) + PIR +
    factor(smoking) + factor(alcohol) + factor(phys_active) + weekday_sleep,
  design = design, family = quasibinomial()
)
m3 <- svyglm(
  dyslipidemia ~ WCUS_cat + age + factor(sex) + factor(race) + factor(education) + PIR +
    factor(smoking) + factor(alcohol) + factor(phys_active) + weekday_sleep +
    BMI + factor(diabetes) + factor(hypertension),
  design = design, family = quasibinomial()
)

# Extract ORs for forest plot
extract_forest_data <- function(model, model_label) {
  cat_names <- c("WCUS_cat0to2", "WCUS_cat2to4", "WCUS_cat4plus")
  cat_labels <- c("0 to <2 h", "2 to <4 h", "≥4 h")
  or <- exp(coef(model))
  ci <- exp(confint(model))
  data.frame(
    model    = model_label,
    category = cat_labels,
    OR       = or[cat_names],
    OR_low   = ci[cat_names, 1],
    OR_high  = ci[cat_names, 2]
  )
}

forest_data <- bind_rows(
  extract_forest_data(m1, "Model 1"),
  extract_forest_data(m2, "Model 2"),
  extract_forest_data(m3, "Model 3")
) %>%
  mutate(
    category = factor(category, levels = c("0 to <2 h", "2 to <4 h", "≥4 h")),
    model    = factor(model, levels = c("Model 3", "Model 2", "Model 1"))
  )

fig3 <- ggplot(forest_data, aes(x = OR, y = model, color = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.25, linewidth = 0.8) +
  geom_point(size = 3) +
  facet_wrap(~category, ncol = 3) +
  scale_x_log10() +
  scale_color_manual(values = c("Model 1" = "#E69F00",
                                "Model 2" = "#56B4E9",
                                "Model 3" = "#009E73")) +
  labs(
    x     = "Odds Ratio (95% CI)",
    y     = NULL,
    color = NULL,
    title = "Association between WCUS Categories and Dyslipidemia"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(hjust = 0.5, size = 11)
  )

ggsave("Figure3_Forest.pdf", fig3, width = 9, height = 4, dpi = 300)
ggsave("Figure3_Forest.png", fig3, width = 9, height = 4, dpi = 300)
cat("Saved: Figure3_Forest.pdf / .png\n")

cat("\nAll figures generated.\n")
sessionInfo()
