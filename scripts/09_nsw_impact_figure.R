#!/usr/bin/env Rscript

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr, tidyr, readr, stringr, ggplot2, scales, 
  haven, lubridate, patchwork, ggpubr
)

testing_path <- "Data/NSW_HCV_notifications_2019_2024_split_age_only.csv"
pbs_path <- "Data/Get data/PBS HCV dataset_March2025.dta"

out_table_csv <- "Results/Impact/Testing_vs_Treatment_AgeCutoffs_2019_2024.csv"
out_plot_cutoffs <- "Results/Impact/Plot_CutoffCoverage_Testing_vs_Treatment_2019_2024.png"
out_plot_breakdown <- "Results/Impact/Plot_AgeBreakdown_Testing_vs_Treatment_2019_2024.png"
out_plot_pyramid <- "Results/Impact/Plot_Pyramid_Testing_vs_Treatment_Age_2019_2024.png"
out_plot_lollipop <- "Results/Impact/Plot_Lollipop_Testing_vs_Treatment_Age_2019_2024.png"
out_plot_overtime <- "Results/Impact/Plot_OverTime_AgeGroupLines_Testing_vs_Treatment_2019_2024.png"
out_plot_stacked <- "Results/Impact/Plot_StackedBars_AgeGroup_Testing_vs_Treatment_2019_2024.png"
out_table_prop_treated <- "Results/Impact/PropTreatedAmongPositive_byAgeGroup_byYear_2019_2024.csv"
out_plot_prop_treated <- "Results/Impact/Plot_PropTreatedAmongPositive_byAgeGroup_byYear_2019_2024.png"
out_table_under_cutoff <- "Results/Impact/Impact_UnderCutoff_NoFibrosisScreen_byYear_2019_2024.csv"
out_plot_under_cutoff <- "Results/Impact/Plot_Impact_UnderCutoff_NoFibrosisScreen_byYear_2019_2024.png"

years_keep <- 2019:2024
cutoffs <- c(35, 40, 45)
pbs_state_code <- 2  # NSW (see pt_st value labels in the PBS dataset)
exclude_age_lower_ge <- 70
fsm_threshold <- 12.5
fsm_threshold_csph <- 20

# --- Helpers ---
parse_age_lower <- function(age_group) {
  age_group <- trimws(age_group)
  lower <- suppressWarnings(as.numeric(str_extract(age_group, "^\\d+")))
  lower
}

make_age_group_levels <- function() {
  c(
    "0-14",
    "15-19",
    "20-24",
    "25-29",
    "30-34",
    "35-39",
    "40-44",
    "45-49",
    "50-54",
    "55-59",
    "60-64",
    "65-69",
    "70-74",
    "75-79",
    "80-84",
    "85+"
  )
}

cut_age_group <- function(age) {
  # Align bins to the NSW CSV labels.
  cut(
    age,
    breaks = c(0, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, Inf),
    right = FALSE,
    labels = make_age_group_levels()
  )
}

# --- Load + reshape NSW testing ---
testing_wide <- read_csv(testing_path, show_col_types = FALSE)

expected_cols <- c("Characteristic")
if (!all(expected_cols %in% names(testing_wide))) {
  stop("Expected column 'Characteristic' not found in testing CSV")
}

# Convert wide year_count/year_percent columns to long.
testing_long <- testing_wide %>%
  pivot_longer(
    cols = -Characteristic,
    names_to = c("year", "metric"),
    names_pattern = "^(\\d{4})_(count|percent)$",
    values_to = "value"
  ) %>%
  mutate(
    year = as.integer(year),
    age_group = as.character(Characteristic)
  ) %>%
  select(year, age_group, metric, value) %>%
  filter(year %in% years_keep) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(
    source = "Testing",
    count = as.numeric(count)
  ) %>%
  select(source, year, age_group, count, percent)

# Exclude minors. The NSW table has 15-19 mixed; to avoid minors entirely we drop that bin.
# This also removes ages 18-19 (adults) due to binning limitations.
testing_adult <- testing_long %>%
  mutate(age_lower = parse_age_lower(age_group)) %>%
  filter(!is.na(age_lower), age_lower >= 20, age_lower < exclude_age_lower_ge)

# Recompute adult-only percent (since provided percent includes minors).
testing_adult <- testing_adult %>%
  group_by(year) %>%
  mutate(
    total_adult = sum(count, na.rm = TRUE),
    adult_percent = if_else(total_adult > 0, count / total_adult, NA_real_)
  ) %>%
  ungroup()

# --- Load PBS treatments and aggregate by year + age group (unique patients) ---
pbs <- read_dta(pbs_path)

needed_pbs_cols <- c("pt_id", "pt_age", "datesupply")
missing_pbs_cols <- setdiff(needed_pbs_cols, names(pbs))
if (length(missing_pbs_cols) > 0) {
  stop("PBS dataset missing required columns: ", paste(missing_pbs_cols, collapse = ", "))
}

pbs_long <- pbs %>%
  mutate(
    year = year(as.Date(datesupply)),
    age = as.numeric(pt_age)
  ) %>%
  filter(pt_st == pbs_state_code) %>%
  filter(year %in% years_keep) %>%
  filter(!is.na(age)) %>%
  mutate(
    age_group = as.character(cut_age_group(age))
  ) %>%
  filter(!is.na(age_group)) %>%
  group_by(year, age_group) %>%
  summarise(count = n_distinct(pt_id), .groups = "drop") %>%
  mutate(source = "Treatment") %>%
  select(source, year, age_group, count)

pbs_adult <- pbs_long %>%
  mutate(age_lower = parse_age_lower(age_group)) %>%
  filter(!is.na(age_lower), age_lower >= 20, age_lower < exclude_age_lower_ge) %>%
  group_by(year) %>%
  mutate(
    total_adult = sum(count, na.rm = TRUE),
    adult_percent = if_else(total_adult > 0, count / total_adult, NA_real_)
  ) %>%
  ungroup()

# --- Cutoff coverage (share captured by age cutoffs) ---
compute_cutoff_coverage_from_bands <- function(df_bands) {
  # df_bands columns: source, year, age_group, count
  df_bands %>%
    mutate(age_lower = parse_age_lower(age_group)) %>%
    filter(!is.na(age_lower)) %>%
    tidyr::expand_grid(cutoff = cutoffs) %>%
    group_by(source, year, cutoff) %>%
    summarise(
      total_adult = sum(count, na.rm = TRUE),
      eligible = sum(count[age_lower >= cutoff], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      pct_eligible = if_else(total_adult > 0, eligible / total_adult, NA_real_),
      cutoff_label = paste0("Age ≥ ", cutoff)
    )
}

cutoff_testing <- compute_cutoff_coverage_from_bands(testing_adult)
cutoff_pbs <- compute_cutoff_coverage_from_bands(pbs_adult)

cutoff_all <- bind_rows(cutoff_testing, cutoff_pbs)

# Wide summary table by source-year
cutoff_table <- cutoff_all %>%
  select(source, year, total_adult, cutoff, eligible, pct_eligible) %>%
  mutate(
    cutoff = paste0("age_ge_", cutoff)
  ) %>%
  pivot_wider(
    names_from = cutoff,
    values_from = c(eligible, pct_eligible)
  ) %>%
  arrange(source, year)

write_csv(cutoff_table, out_table_csv)

# --- Plot 1: cutoff coverage comparison ---
plot_cutoffs <- cutoff_all %>%
  mutate(
    cutoff_label = factor(cutoff_label, levels = paste0("Age ≥ ", cutoffs)),
    year = as.integer(year)
  ) %>%
  ggplot(aes(x = year, y = pct_eligible, color = source, group = source)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~ cutoff_label, ncol = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = years_keep) +
  labs(
    title = "Impact of age cutoffs: notifications vs treatment (NSW, 2019–2024)",
    subtitle = "Adult-only (excludes <20 due to 15–19 bin containing minors)",
    x = "Year",
    y = "Proportion captured by cutoff",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(out_plot_cutoffs, plot_cutoffs, width = 9, height = 6, dpi = 300)

# --- Plot 2: age breakdown over time (testing vs treatment) ---
adult_age_levels_u70 <- setdiff(
  make_age_group_levels(),
  c("0-14", "15-19", "70-74", "75-79", "80-84", "85+")
)

breakdown <- bind_rows(
  testing_adult %>% select(source, year, age_group, count),
  pbs_adult %>% select(source, year, age_group, count)
) %>%
  mutate(
    age_group = factor(age_group, levels = adult_age_levels_u70),
    year = as.integer(year)
  )

plot_breakdown <- breakdown %>%
  ggplot(aes(x = year, y = count, fill = age_group)) +
  geom_col(position = "fill", color = NA) +
  facet_wrap(~ source, ncol = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = years_keep) +
  labs(
    title = "Age distribution over time: notifications vs treatment (NSW, 2019–2024)",
    subtitle = "Proportions within each year (adult-only; excludes <20)",
    x = "Year",
    y = "Within-year proportion",
    fill = "Age group"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

ggsave(out_plot_breakdown, plot_breakdown, width = 10, height = 8, dpi = 300)

# --- Plot 3: population pyramid (Testing vs Treatment; years combined; adult-only) ---
pyramid <- bind_rows(
  testing_adult %>% select(source, age_group, count),
  pbs_adult %>% select(source, age_group, count)
) %>%
  group_by(source, age_group) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(source) %>%
  mutate(
    total = sum(count, na.rm = TRUE),
    pct = if_else(total > 0, count / total, NA_real_)
  ) %>%
  ungroup() %>%
  mutate(
    # Reverse so youngest appears at the top of the y-axis.
    age_group = factor(age_group, levels = rev(adult_age_levels_u70)),
    x = if_else(source == "Testing", -pct, pct),
    source = factor(source, levels = c("Testing", "Treatment")),
    pct_label = percent(pct, accuracy = 0.1),
    label_x = if_else(source == "Testing", x - 0.01, x + 0.01)
  )

plot_pyramid <- pyramid %>%
  ggplot(aes(y = age_group, x = x, fill = source)) +
  geom_col(width = 0.9, alpha = 0.9) +
  geom_text(
    aes(x = label_x, label = pct_label, hjust = if_else(source == "Testing", 1, 0)),
    size = 3,
    color = "black"
  ) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.6) +
  scale_x_continuous(
    labels = function(z) percent(abs(z), accuracy = 1),
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = 0.10)
  ) +
  labs(
    title = "Population pyramid: notifications vs treatment age distribution (NSW, 2019–2024)",
    subtitle = "Adult-only (excludes <20 because the testing table has a 15–19 bin containing minors)",
    x = "Percent (Testing left, Treatment right)",
    y = "Age group",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(out_plot_pyramid, plot_pyramid, width = 10, height = 7, dpi = 300)

# --- Plot 4: lollipop contrast (Testing vs Treatment; years combined; adult-only) ---
lollipop_long <- pyramid %>%
  select(source, age_group, pct)

lollipop_wide <- lollipop_long %>%
  pivot_wider(names_from = source, values_from = pct)

lollipop_points <- lollipop_wide %>%
  pivot_longer(cols = c(Testing, Treatment), names_to = "source", values_to = "pct") %>%
  mutate(source = factor(source, levels = c("Testing", "Treatment")))

plot_lollipop <- ggplot() +
  geom_segment(
    data = lollipop_wide,
    aes(y = age_group, yend = age_group, x = Testing, xend = Treatment),
    linewidth = 0.8,
    color = "grey50"
  ) +
  geom_point(
    data = lollipop_points,
    aes(x = pct, y = age_group, color = source),
    size = 2.6
  ) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "Notifications vs treatment by age group (NSW, 2019–2024)",
    subtitle = "Lollipop plot: within-source percentages; adult-only (excludes <20)",
    x = "Percent",
    y = "Age group",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(out_plot_lollipop, plot_lollipop, width = 9, height = 7, dpi = 300)

# --- Plot 5: over-time lines by age group (Testing vs Treatment panels) ---
overtime <- bind_rows(
  testing_adult %>% select(source, year, age_group, count),
  pbs_adult %>% select(source, year, age_group, count)
) %>%
  mutate(
    year = as.integer(year),
    age_group = as.character(age_group)
  ) %>%
  filter(age_group %in% adult_age_levels_u70) %>%
  tidyr::complete(
    source,
    year = years_keep,
    age_group = adult_age_levels_u70,
    fill = list(count = 0)
  ) %>%
  group_by(source, year) %>%
  mutate(
    total = sum(count, na.rm = TRUE),
    pct = if_else(total > 0, count / total, NA_real_)
  ) %>%
  ungroup() %>%
  mutate(
    age_group = factor(age_group, levels = rev(adult_age_levels_u70)),
    source = factor(source, levels = c("Testing", "Treatment"))
  )

plot_overtime <- overtime %>%
  ggplot(aes(x = year, y = pct, color = age_group, group = age_group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  facet_wrap(~ source, ncol = 2) +
  scale_x_continuous(breaks = years_keep) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "Age distribution over time: notifications vs treatment (NSW, 2019–2024)",
    subtitle = "Years on x-axis; within-year percent on y-axis; adult-only (excludes <20)",
    x = "Year",
    y = "Within-year percent",
    color = "Age group"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

ggsave(out_plot_overtime, plot_overtime, width = 12, height = 7, dpi = 300)

# --- Plot 6: stacked bars (within-year %) by age group ---
plot_stacked <- overtime %>%
  ggplot(aes(x = year, y = pct, fill = age_group)) +
  geom_col(width = 0.85) +
  facet_wrap(~ source, ncol = 2) +
  scale_x_continuous(breaks = years_keep) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Age distribution over time: stacked bars (NSW, 2019–2024)",
    subtitle = "Within-year percent by age group; adult-only (excludes <20)",
    x = "Year",
    y = "Within-year percent",
    fill = "Age group"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

ggsave(out_plot_stacked, plot_stacked, width = 12, height = 7, dpi = 300)

# --- Treated among positives (by year x age group) ---
testing_pos <- testing_adult %>%
  transmute(year = as.integer(year), age_group = as.character(age_group), positives = as.numeric(count))

treat_nsw <- pbs_adult %>%
  transmute(year = as.integer(year), age_group = as.character(age_group), treated = as.numeric(count))

prop_treated <- full_join(testing_pos, treat_nsw, by = c("year", "age_group")) %>%
  mutate(
    positives = coalesce(positives, 0),
    treated = coalesce(treated, 0),
    treated_over_positives = if_else(positives > 0, treated / positives, NA_real_),
    age_group = factor(age_group, levels = rev(adult_age_levels_u70))
  ) %>%
  filter(year %in% years_keep, !is.na(age_group)) %>%
  arrange(year, age_group)

write_csv(
  prop_treated %>%
    mutate(
      age_group = as.character(age_group)
    ) %>%
    select(year, age_group, positives, treated, treated_over_positives),
  out_table_prop_treated
)

plot_prop_treated <- prop_treated %>%
  ggplot(aes(x = year, y = age_group, fill = treated_over_positives)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(
    aes(label = if_else(is.na(treated_over_positives), "", percent(treated_over_positives, accuracy = 1))),
    size = 3
  ) +
  scale_x_continuous(breaks = years_keep) +
  scale_fill_gradient(
    low = "white",
    high = "steelblue",
    labels = percent_format(accuracy = 1),
    na.value = "grey90"
  ) +
  labs(
    title = "Treated / positive notifications (NSW, 2019–2024)",
    subtitle = "PBS unique patients treated divided by positive notifications; can exceed 100% (adult-only; excludes <20)",
    x = "Year",
    y = "Age group",
    fill = "Treated / positive"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

ggsave(out_plot_prop_treated, plot_prop_treated, width = 11, height = 7, dpi = 300)

message("Wrote: ", out_table_csv)
message("Wrote: ", out_plot_cutoffs)
message("Wrote: ", out_plot_breakdown)
message("Wrote: ", out_plot_pyramid)
message("Wrote: ", out_plot_lollipop)
message("Wrote: ", out_plot_overtime)
message("Wrote: ", out_plot_stacked)
message("Wrote: ", out_table_prop_treated)
message("Wrote: ", out_plot_prop_treated)

# --- Impact of age cutoffs (reverse): under-cutoff would not get fibrosis screening ---
# Use NSW positive-notification counts as the population needing fibrosis assessment.
# Estimate the number with advanced fibrosis (FSM >= threshold) by age band using merged cohort.
cohort_fsm_by_ageband <- read_csv("Data/merged_data.csv", show_col_types = FALSE) %>%
  filter(!is.na(age), !is.na(fsm)) %>%
  mutate(
    age_group = as.character(cut_age_group(as.numeric(age))),
    age_lower = parse_age_lower(age_group),
    fsm_high = as.numeric(fsm >= fsm_threshold),
    fsm_csph = as.numeric(fsm > fsm_threshold_csph)
  ) %>%
  filter(!is.na(age_group), !is.na(age_lower), age_lower >= 20, age_lower < exclude_age_lower_ge) %>%
  group_by(age_group) %>%
  summarise(
    cohort_n = n(),
    p_fsm_high = mean(fsm_high, na.rm = TRUE),
    p_fsm_csph = mean(fsm_csph, na.rm = TRUE),
    .groups = "drop"
  )

testing_pos_by_year_age <- testing_adult %>%
  transmute(
    year = as.integer(year),
    age_group = as.character(age_group),
    age_lower = parse_age_lower(age_group),
    positives = as.numeric(count)
  ) %>%
  filter(!is.na(age_lower), age_lower >= 20, age_lower < exclude_age_lower_ge)

under_cutoff <- testing_pos_by_year_age %>%
  left_join(cohort_fsm_by_ageband, by = "age_group") %>%
  mutate(
    p_fsm_high = coalesce(p_fsm_high, 0),
    p_fsm_csph = coalesce(p_fsm_csph, 0),
    expected_fsm_high = positives * p_fsm_high,
    expected_fsm_csph = positives * p_fsm_csph
  ) %>%
  group_by(year) %>%
  mutate(total_positives = sum(positives, na.rm = TRUE)) %>%
  ungroup() %>%
  tidyr::expand_grid(cutoff = cutoffs) %>%
  group_by(year, cutoff) %>%
  summarise(
    total_positives = first(total_positives),
    under_cutoff_positives = sum(positives[age_lower < cutoff], na.rm = TRUE),
    pct_under_cutoff = if_else(total_positives > 0, under_cutoff_positives / total_positives, NA_real_),
    expected_missed_fsm_high = sum(expected_fsm_high[age_lower < cutoff], na.rm = TRUE),
    expected_missed_fsm_csph = sum(expected_fsm_csph[age_lower < cutoff], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cutoff_label = factor(paste0("Age < ", cutoff), levels = paste0("Age < ", cutoffs))
  ) %>%
  filter(year %in% years_keep)

write_csv(under_cutoff, out_table_under_cutoff)

plot_bars <- under_cutoff %>%
  ggplot(aes(x = year, y = pct_under_cutoff)) +
  geom_col(aes(fill = cutoff), width = 0.85) +
  geom_text(
    aes(label = percent(pct_under_cutoff, accuracy = 1)),
    vjust = -0.5,
    size = 4
  ) +
  facet_wrap(~ cutoff_label, ncol = 3) +
  scale_x_continuous(breaks = years_keep) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.10))) +
  labs(
    title = "Under-cutoff excluded from fibrosis screening (NSW notifications, 2019–2024)",
    x = "Year",
    y = "Percent avoiding fibrosis assessment"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none"
  )

table_data_wide <- under_cutoff %>%
  select(year, cutoff_label, under_cutoff_positives, expected_missed_fsm_high) %>%
  pivot_longer(c(under_cutoff_positives, expected_missed_fsm_high), names_to = "Metric", values_to = "value") %>%
  mutate(
    Metric = if_else(Metric == "under_cutoff_positives", "Not screened", "Missed cirrhosis"),
    Metric = factor(Metric, levels = c("Not screened", "Missed cirrhosis")),
    value = format(round(value), big.mark = ",")
  ) %>%
  pivot_wider(names_from = "year", values_from = "value") %>%
  arrange(cutoff_label, Metric)

tbls <- lapply(levels(under_cutoff$cutoff_label), function(cut) {
  dt <- table_data_wide %>%
    filter(cutoff_label == cut) %>%
    select(-cutoff_label) %>%
    rename(` ` = Metric)
  
  # Remove year column names so they don't print in the table (years are in the main plot directly above)
  colnames(dt)[-1] <- rep("", ncol(dt) - 1)
  
  # Return a ggtexttable object
  ggtexttable(
    dt, 
    rows = NULL, 
    theme = ttheme("minimal", base_size = 10, padding = unit(c(2, 4), "mm"))
  )
})

plot_table <- wrap_plots(tbls, ncol = 3)

plot_under_cutoff <- plot_bars / plot_table + plot_layout(heights = c(4, 1))

plot_under_cutoff

ggsave(out_plot_under_cutoff, plot_under_cutoff, width = 12, height = 7, dpi = 300)

message("Wrote: ", out_table_under_cutoff)
message("Wrote: ", out_plot_under_cutoff)

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)

under_cutoff <- under_cutoff %>%
  mutate(
    cutoff_label = factor(cutoff_label, levels = c("Age < 35", "Age < 40", "Age < 45")),
    year = factor(year, levels = years_keep)
  )

under_cutoff <- under_cutoff %>%
  mutate(
    cutoff_label = case_match(
      cutoff_label,
      "Age < 35" ~ "Age > 35 Strategy",
      "Age < 40" ~ "Age > 40 Strategy",
      "Age < 45" ~ "Age > 45 Strategy"
    )
  )

# -------------------------
# Top bar plot
# -------------------------
plot_bars <- ggplot(under_cutoff, aes(x = year, y = pct_under_cutoff)) +
  geom_col(width = 0.85, fill = "#3B6EA5") +
  geom_text(
    aes(label = percent(pct_under_cutoff, accuracy = 1)),
    vjust = -0.4,
    size = 4
  ) +
  facet_grid(~ cutoff_label) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    title = "Under-cutoff excluded from fibrosis screening (NSW notifications, 2019–2024)",
    x = "Year",
    y = "Percent avoiding fibrosis assessment"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 5.5, 0, 5.5)
  )

# -------------------------
# Bottom table data
# Top row in table = Not screened
# Bottom row in table = Missed cirrhosis
# -------------------------
table_df <- under_cutoff %>%
  transmute(
    cutoff_label,
    year,
    `Not screened` = format(round(under_cutoff_positives), big.mark = ","),
    `Missed cirrhosis` = format(round(expected_missed_fsm_high), big.mark = ","),
    `Missed CSPH` = format(round(expected_missed_fsm_csph), big.mark = ",")
  ) %>%
  pivot_longer(
    cols = c(`Not screened`, `Missed cirrhosis`, `Missed CSPH`),
    names_to = "metric",
    values_to = "label"
  ) %>%
  mutate(
    metric = factor(metric, levels = c("Not screened", "Missed cirrhosis", "Missed CSPH")),
    y = case_match(
      metric,
      "Not screened" ~ 3,
      "Missed cirrhosis" ~ 2,
      "Missed CSPH" ~ 1
    )
  )

# -------------------------
# Bottom table plot with borders
# Years and x-axis title shown on TOP of this panel
# -------------------------
plot_table <- ggplot(table_df, aes(x = year, y = y)) +
  geom_tile(
    width = 1,
    height = 1,
    fill = "white",
    color = "grey70",
    linewidth = 0.4
  ) +
  geom_text(aes(label = label), size = 3.6) +
  facet_grid(~ cutoff_label) +
  scale_y_continuous(
    limits = c(0.5, 3.5),
    breaks = c(1, 2, 3),
    labels = c("", "", "")
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.title.x.top = element_text(size = 11, margin = margin(b = 8)),
    
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.title.x.bottom = element_blank(),
    
    plot.margin = margin(0, 5.5, 5.5, 5.5)
  ) +
  scale_x_discrete(position = "top")

# -------------------------
# Combine
# -------------------------
plot_under_cutoff <- plot_bars / plot_table +
  plot_layout(heights = c(8, 2.0))

plot_under_cutoff

ggsave("Results/Impact/Plot_Impact_UnderCutoff_NoFibrosisScreen_byYear_2019_2024.pdf", plot_under_cutoff, width = 12, height = 7, dpi = 300)
ggsave(out_plot_under_cutoff, plot_under_cutoff, width = 12, height = 7, dpi = 300)

