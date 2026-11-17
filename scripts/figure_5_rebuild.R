#!/usr/bin/env Rscript

# Run from the repository root after placing the private inputs in Data/.
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tidyr, readr, stringr, ggplot2, scales, patchwork)

cohort_path <- "Data/merged_data.csv"
notifications_path <- "Data/NSW_HCV_notifications_2019_2024_split_age_only.csv"
output_dir <- "outputs/figures"
stopifnot(file.exists(cohort_path), file.exists(notifications_path))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

years <- 2019:2024
cutoffs <- c(35, 40, 45)

# Derive five-year prevalences from the included participant cohort.
cohort <- read_csv(cohort_path, show_col_types = FALSE) %>%
  filter(!is.na(age), !is.na(fsm)) %>%
  mutate(
    age_lower = 5 * floor(age / 5),
    age_group = paste0(age_lower, "-", age_lower + 4),
    cirrhosis = as.numeric(fsm >= 12.5),
    csph = as.numeric(fsm > 20)
  ) %>%
  filter(age_lower >= 20, age_lower < 70)

prevalence <- cohort %>%
  group_by(age_group) %>%
  summarise(
    p_cirrhosis = mean(cirrhosis),
    p_csph = mean(csph),
    .groups = "drop"
  )

# The NSW file supplies five-year age bands and annual notification counts.
notifications <- read_csv(notifications_path, show_col_types = FALSE) %>%
  pivot_longer(
    cols = -Characteristic,
    names_to = c("year", "metric"),
    names_pattern = "^(\\d{4})_(count|percent)$",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  transmute(
    year = as.integer(year),
    age_group = Characteristic,
    age_lower = as.numeric(str_extract(Characteristic, "^\\d+")),
    count = as.numeric(count)
  ) %>%
  filter(year %in% years, age_lower >= 20, age_lower < 70) %>%
  left_join(prevalence, by = "age_group")

stopifnot(!anyNA(notifications[, c("count", "p_cirrhosis", "p_csph")]))

# Estimate cases whose baseline identification would be deferred.
impact <- notifications %>%
  mutate(
    expected_cirrhosis = count * p_cirrhosis,
    expected_csph = count * p_csph
  ) %>%
  crossing(cutoff = cutoffs) %>%
  group_by(year, cutoff) %>%
  summarise(
    total = sum(count),
    not_screened = sum(count[age_lower < cutoff]),
    missed_cirrhosis = sum(expected_cirrhosis[age_lower < cutoff]),
    missed_csph = sum(expected_csph[age_lower < cutoff]),
    .groups = "drop"
  ) %>%
  mutate(
    pct_not_screened = not_screened / total,
    cutoff_label = factor(paste0("Age < ", cutoff), levels = paste0("Age < ", cutoffs))
  )

stopifnot(sum(impact$not_screened[impact$cutoff == 35]) == 2826)
stopifnot(sum(impact$not_screened[impact$cutoff == 40]) == 3709)
stopifnot(sum(impact$not_screened[impact$cutoff == 45]) == 4523)
stopifnot(identical(as.integer(round(tapply(impact$missed_csph, impact$cutoff, sum))), c(19L, 45L, 91L)))

# Reserve identical positions for row labels and years in both plot sections.
year_levels <- as.character(years)
slots <- c("label_1", "label_2", year_levels)
impact <- impact %>% mutate(year = factor(as.character(year), levels = slots))

bars <- ggplot(impact, aes(year, pct_not_screened)) +
  geom_col(width = 0.85, fill = "#3B6EA5") +
  geom_text(aes(label = percent(pct_not_screened, accuracy = 1)), vjust = -0.5, size = 4) +
  facet_grid(~ cutoff_label) +
  scale_x_discrete(limits = slots, labels = c("", "", year_levels), drop = FALSE,
                   expand = expansion(add = 0.5)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.10))) +
  labs(x = "Year", y = "Percent avoiding fibrosis assessment") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 5.5, 0, 5.5))

table_values <- impact %>%
  transmute(
    cutoff_label, year,
    `Not screened` = format(round(not_screened), big.mark = ","),
    `Missed cirrhosis` = format(round(missed_cirrhosis), big.mark = ","),
    `Missed CSPH` = format(round(missed_csph), big.mark = ",")
  ) %>%
  pivot_longer(c(`Not screened`, `Missed cirrhosis`, `Missed CSPH`),
               names_to = "metric", values_to = "label") %>%
  mutate(metric = factor(metric,
                         levels = c("Not screened", "Missed cirrhosis", "Missed CSPH")),
         row = match(metric, rev(levels(metric))))

row_labels <- expand_grid(
  cutoff_label = levels(impact$cutoff_label),
  metric = levels(table_values$metric)
) %>%
  mutate(row = match(metric, rev(levels(table_values$metric))),
         year = factor("label_2", levels = slots))

table_plot <- ggplot(table_values, aes(year, row)) +
  geom_tile(width = 1, height = 1, fill = "white", color = "grey70", linewidth = 0.4) +
  geom_text(aes(label = label), size = 3.2) +
  geom_text(data = row_labels, aes(year, row, label = metric), inherit.aes = FALSE,
            hjust = 1, size = 2.5) +
  facet_grid(~ cutoff_label) +
  scale_x_discrete(limits = slots, labels = rep("", length(slots)), drop = FALSE,
                   expand = expansion(add = 0.5)) +
  scale_y_continuous(limits = c(0.5, 3.5), breaks = 1:3, labels = c("", "", "")) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), plot.margin = margin(0, 5.5, 5.5, 5.5))

figure <- bars / table_plot + plot_layout(heights = c(8, 2))
ggsave(file.path(output_dir, "Figure_5_600dpi.eps"), figure,
       width = 12, height = 7, device = "eps")
ggsave(file.path(output_dir, "Figure_5_600dpi.png"), figure,
       width = 12, height = 7, dpi = 600)
