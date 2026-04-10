# Asssociation between weekend catch-up sleep and dyslipidemia among U.S. adults

## Overview

R analysis code for the manuscript:

> **"Association between weekend catch-up sleep and dyslipidemia among U.S. adults: A cross-sectional study of NHANES 2017–March 2020"**  
> Jie He  
> *Submitted to PLOS ONE*

## Data Source

All data are publicly available from the U.S. Centers for Disease Control and Prevention (CDC) National Health and Nutrition Examination Survey (NHANES):

🔗 https://www.cdc.gov/nchs/nhanes/

Specifically, the 2017–March 2020 Pre-Pandemic cycle files (prefix `P_`) were used.

## R Version and Packages

- R version 4.5.3 (2026-03-11)
- `tidyverse` — data manipulation
- `haven` — reading NHANES .xpt files
- `survey` — complex survey-weighted analysis
- `rms` — restricted cubic splines
- `ggplot2` — figure generation

## Scripts

| Script | Description |
|--------|-------------|
| `01_data_preparation.R` | Load and merge NHANES XPT files, define variables, apply exclusion criteria, save analytic dataset |
| `02_main_analysis.R` | Survey-weighted logistic regression (Models 1–3), stratified analyses, individual lipid components |
| `03_sensitivity_analyses.R` | Sensitivity analyses S1–S4 |
| `04_figures.R` | Generate Figure 2 (RCS dose-response) and Figure 3 (forest plot) |

## How to Reproduce

1. Download the following NHANES 2017–March 2020 XPT files from https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?Cycle=2017-2020:
   - P_DEMO.xpt, P_SLQ.xpt, P_TCHOL.xpt, P_TRIGLY.xpt, P_HDL.xpt
   - P_BMX.xpt, P_GHB.xpt, P_GLU.xpt, P_BPQ.xpt, P_BPXO.xpt
   - P_DIQ.xpt, P_SMQ.xpt, P_ALQ.xpt, P_PAQ.xpt

2. Place all XPT files in your working directory.

3. Run scripts in order:
```r
source("01_data_preparation.R")   # Creates nhanes_analytic_dataset.csv
source("02_main_analysis.R")      # Main results
source("03_sensitivity_analyses.R")
source("04_figures.R")
```

## Study Design

- **Design:** Cross-sectional
- **Data:** NHANES 2017–March 2020 Pre-Pandemic
- **Sample:** U.S. adults ≥18 years, fasting subsample (N = 2,979)
- **Exposure:** Weekend catch-up sleep (WCUS = weekend − weekday sleep duration)
- **Outcome:** Dyslipidemia (ATP III / AHA criteria)
- **Analysis:** Survey-weighted logistic regression using the `survey` package with Taylor series linearization

## License

This code is released under the [MIT License](LICENSE).

## Contact

For questions about the code, please open a GitHub Issue or contact the corresponding author, Jie He, at hejiewljyb@163.com
