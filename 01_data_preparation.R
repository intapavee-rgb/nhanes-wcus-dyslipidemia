# =============================================================================
# Script 01: Data Preparation
# Project: Association between weekend catch-up sleep and dyslipidemia among U.S. adults: A cross-sectional study of NHANES 2017–March 2020
# Data: NHANES 2017–March 2020 Pre-Pandemic
# Author: Paveena Intaraksa
# =============================================================================

# Install packages if needed:
# install.packages(c("tidyverse", "haven", "survey"))

library(tidyverse)
library(haven)
library(survey)

# =============================================================================
# SECTION 1: LOAD RAW NHANES XPT FILES
# =============================================================================
# Download files from: https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?Cycle=2017-2020
# Place all .xpt files in your working directory

# Demographics
demo  <- read_xpt("P_DEMO.xpt") %>%
  select(SEQN, RIDSTATR, RIAGENDR, RIDAGEYR, RIDRETH3, DMDEDUC2,
         RIDEXPRG, INDFMPIR, WTMECPRP, WTSAFPRP, SDMVPSU, SDMVSTRA)

# Sleep
slq   <- read_xpt("P_SLQ.xpt")  %>% select(SEQN, SLD012, SLD013)

# Lipids
tchol  <- read_xpt("P_TCHOL.xpt")  %>% select(SEQN, LBXTC)
trigly <- read_xpt("P_TRIGLY.xpt") %>% select(SEQN, LBXTR, LBDLDL)
hdl    <- read_xpt("P_HDL.xpt")    %>% select(SEQN, LBDHDD)

# Clinical
bmx  <- read_xpt("P_BMX.xpt")  %>% select(SEQN, BMXBMI)
ghb  <- read_xpt("P_GHB.xpt")  %>% select(SEQN, LBXGH)
glu  <- read_xpt("P_GLU.xpt")  %>% select(SEQN, LBXGLU)
bpq  <- read_xpt("P_BPQ.xpt")  %>% select(SEQN, BPQ040A, BPQ100D)
bpx  <- read_xpt("P_BPXO.xpt") %>% select(SEQN, BPXOSY1, BPXODI1, BPXOSY2, BPXODI2)
diq  <- read_xpt("P_DIQ.xpt")  %>% select(SEQN, DIQ010, DIQ070)

# Lifestyle
smq  <- read_xpt("P_SMQ.xpt")  %>% select(SEQN, SMQ020, SMQ040)
alq  <- read_xpt("P_ALQ.xpt")  %>% select(SEQN, ALQ111, ALQ121)
paq  <- read_xpt("P_PAQ.xpt")  %>% select(SEQN, PAQ610, PAD615, PAQ655, PAD660)

# =============================================================================
# SECTION 2: MERGE ALL FILES
# =============================================================================

nhanes_raw <- demo %>%
  left_join(slq,    by = "SEQN") %>%
  left_join(tchol,  by = "SEQN") %>%
  left_join(trigly, by = "SEQN") %>%
  left_join(hdl,    by = "SEQN") %>%
  left_join(bmx,    by = "SEQN") %>%
  left_join(ghb,    by = "SEQN") %>%
  left_join(glu,    by = "SEQN") %>%
  left_join(bpq,    by = "SEQN") %>%
  left_join(bpx,    by = "SEQN") %>%
  left_join(diq,    by = "SEQN") %>%
  left_join(smq,    by = "SEQN") %>%
  left_join(alq,    by = "SEQN") %>%
  left_join(paq,    by = "SEQN")

cat("Raw merged dataset:", nrow(nhanes_raw), "rows\n")

# =============================================================================
# SECTION 3: VARIABLE DEFINITIONS AND RECODING
# =============================================================================

nhanes <- nhanes_raw %>%

  # --- Eligibility filters ---
  filter(RIDSTATR == 2) %>%                         # Examined participants only
  filter(RIDAGEYR >= 18) %>%                        # Adults 18+
  filter(is.na(RIDEXPRG) | RIDEXPRG != 1) %>%      # Exclude pregnant women

  mutate(

    # --- Weekend Catch-Up Sleep (WCUS) ---
    # SLD012 = self-reported weekday sleep duration (hours)
    # SLD013 = self-reported weekend sleep duration (hours)
    weekday_sleep = SLD012,
    weekend_sleep = SLD013,
    WCUS = weekend_sleep - weekday_sleep,

    # Categorical WCUS (primary analysis)
    WCUS_cat = factor(case_when(
      WCUS <  0              ~ "<0",
      WCUS >= 0 & WCUS <  2 ~ "0to2",
      WCUS >= 2 & WCUS <  4 ~ "2to4",
      WCUS >= 4              ~ "4plus"
    ), levels = c("<0", "0to2", "2to4", "4plus")),

    # Ordinal WCUS for P-for-trend
    WCUS_ordinal = as.numeric(WCUS_cat),

    # --- Dyslipidemia (ATP III / AHA criteria) ---
    high_TC  = (!is.na(LBXTC))  & LBXTC  >= 240,
    high_LDL = (!is.na(LBDLDL)) & LBDLDL >= 160,
    low_HDL  = case_when(
      RIAGENDR == 1 ~ (!is.na(LBDHDD)) & LBDHDD < 40,   # Males
      RIAGENDR == 2 ~ (!is.na(LBDHDD)) & LBDHDD < 50    # Females
    ),
    high_TG   = (!is.na(LBXTR))  & LBXTR  >= 200,
    lipid_med = (!is.na(BPQ100D)) & BPQ100D == 1,
    dyslipidemia = as.integer(high_TC | high_LDL | low_HDL | high_TG | lipid_med),

    # --- Sociodemographic covariates ---
    age = RIDAGEYR,
    sex = RIAGENDR,   # 1 = Male, 2 = Female

    race = case_when(
      RIDRETH3 == 3 ~ 1,   # Non-Hispanic White
      RIDRETH3 == 4 ~ 2,   # Non-Hispanic Black
      RIDRETH3 %in% c(1,2) ~ 3,  # Hispanic
      RIDRETH3 == 6 ~ 4,   # Non-Hispanic Asian
      RIDRETH3 == 7 ~ 5    # Other/Multiracial
    ),

    education = case_when(
      DMDEDUC2 %in% c(1, 2) ~ 1,   # Less than high school
      DMDEDUC2 == 3          ~ 2,   # High school graduate / GED
      DMDEDUC2 == 4          ~ 3,   # Some college / associate
      DMDEDUC2 == 5          ~ 4    # College graduate or above
    ),

    PIR = INDFMPIR,   # Poverty-income ratio (continuous)

    # --- Lifestyle covariates ---
    smoking = case_when(
      SMQ020 == 2                      ~ 0,   # Never
      SMQ020 == 1 & SMQ040 == 3        ~ 1,   # Former
      SMQ020 == 1 & SMQ040 %in% c(1,2) ~ 2   # Current
    ),

    alcohol = case_when(
      ALQ111 == 2 | is.na(ALQ111) ~ 0,                  # Non-drinker
      ALQ111 == 1 & ALQ121 <= 52  ~ 1,                  # Moderate
      ALQ111 == 1 & ALQ121 > 52   ~ 2                   # Heavy
    ),

    # Physical activity: MET-min/week (moderate + vigorous)
    pa_mod  = ifelse(!is.na(PAQ610) & !is.na(PAD615), PAQ610 * PAD615, 0),
    pa_vig  = ifelse(!is.na(PAQ655) & !is.na(PAD660), PAQ655 * PAD660 * 2, 0),
    pa_total = pa_mod + pa_vig,
    phys_active = as.integer(pa_total >= 150),  # 1 = meets guidelines

    # --- Clinical covariates ---
    BMI = BMXBMI,

    diabetes = as.integer(
      (!is.na(DIQ010) & DIQ010 == 1) |
      (!is.na(LBXGH)  & LBXGH  >= 6.5) |
      (!is.na(LBXGLU) & LBXGLU >= 126) |
      (!is.na(DIQ070) & DIQ070 == 1)
    ),

    # Average systolic and diastolic BP across two readings
    SBP = rowMeans(cbind(BPXOSY1, BPXOSY2), na.rm = TRUE),
    DBP = rowMeans(cbind(BPXODI1, BPXODI2), na.rm = TRUE),

    hypertension = as.integer(
      (!is.na(SBP) & SBP >= 130) |
      (!is.na(DBP) & DBP >= 80)  |
      (!is.na(BPQ040A) & BPQ040A == 1)
    )
  )

# =============================================================================
# SECTION 4: EXCLUSIONS AND COMPLETE-CASE ANALYSIS
# =============================================================================

cat("Before exclusions:", nrow(nhanes), "\n")

# Require non-missing values on all key variables
nhanes_complete <- nhanes %>%
  filter(
    !is.na(WCUS_cat),
    !is.na(dyslipidemia),
    !is.na(age),
    !is.na(sex),
    !is.na(race),
    !is.na(education),
    !is.na(PIR),
    !is.na(smoking),
    !is.na(alcohol),
    !is.na(phys_active),
    !is.na(weekday_sleep),
    !is.na(BMI),
    !is.na(diabetes),
    !is.na(hypertension),
    !is.na(WTSAFPRP), WTSAFPRP > 0   # Must have valid fasting weight
  )

cat("Final analytic sample:", nrow(nhanes_complete), "\n")

# =============================================================================
# SECTION 5: SAVE ANALYTIC DATASET
# =============================================================================

write.csv(nhanes_complete, "nhanes_analytic_dataset.csv", row.names = FALSE)
cat("Saved: nhanes_analytic_dataset.csv\n")

# Session info for reproducibility
sessionInfo()
