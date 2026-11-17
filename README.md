# Treat First analysis

R scripts and selected final figures for the individual participant data analysis of age-guided fibrosis assessment before hepatitis C treatment. The study data are governed by contributing-study agreements and are not in this repository.

## Reproduce locally

Use R 4.5 or later. Run commands from the repository root. The scripts use `pacman::p_load()` to install and load their R packages. Install Pandoc if rendering the R Markdown reports.

Place authorized private inputs under `Data/` locally. This directory is ignored by Git. The core scripts expect the study CSV files used by `scripts/01_merge_data.Rmd`; downstream scripts use `Data/merged_data_all.csv`, `Data/merged_data.csv`, and `Data/study_data.csv`. The NSW Figure 5 script also needs `Data/NSW_HCV_notifications_2019_2024_split_age_only.csv`. The risk-of-bias script needs `Data/Risk of bias.xlsx`; PBS analyses need the authorized PBS `.dta` file under `Data/Get data/`.

Run the numbered scripts in order. `01_merge_data.Rmd` creates the harmonized analysis files in the private `Data/` directory. The remaining scripts create results under `Results/`, which is also ignored by Git. `scripts/figure_5_rebuild.R` reproduces the public Figure 5 EPS and PNG from the local merged cohort and NSW notification counts:

```sh
Rscript scripts/figure_5_rebuild.R
```

The Figure 5 calculation defines clinically significant portal hypertension as liver stiffness **greater than 20 kPa**. Its published bars are the share of NSW adult notifications below ages 35, 40, and 45; the bottom rows show estimated counts. The EPS outputs are vector figures and the PNG is exported at 600 dpi.

## Repository contents

- `scripts/`: the analysis notebooks and the Figure 5 rebuild script.
- `outputs/figures/`: five manuscript EPS figures without an overall title, plus a Figure 5 PNG preview.

No raw participant records, harmonized datasets, source workbooks, or CSV result tables are tracked. The figures are aggregate outputs. Check the manuscript and supplement against newly generated outputs before submission.
