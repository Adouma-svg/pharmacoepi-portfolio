# ==============================================================================
# Study 1: Real-World Treatment Patterns in Advanced NSCLC
# Script:   01_analysis_treatment_patterns.R
# Purpose:  Describe treatment patterns, time-to-treatment, and sequencing
#           using methods standard in RWE/pharmacoepidemiology
# Concepts: Table 1, Sankey-style flow, time-to-treatment, logistic regression
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tableone)
  library(gtsummary)
  library(ggalluvial)
  library(scales)
  library(forcats)
  library(stringr)
  library(patchwork)
})

# ── 0. Load data ──────────────────────────────────────────────────────────────
df <- read.csv("study1_treatment_patterns/data/nsclc_cohort.csv",
               stringsAsFactors = FALSE) %>%
  mutate(
    index_date    = as.Date(index_date),
    start_date_1L = as.Date(start_date_1L),
    # Create analysis-ready factor variables
    age_cat   = cut(age, c(29,54,64,74,90), labels = c("<55","55-64","65-74","75+")),
    stage_grp = ifelse(stage %in% c("IIIA","IIIB"), "Stage III", "Stage IV"),
    year_dx   = as.integer(format(index_date, "%Y")),
    pdl1_known = ifelse(pdl1_expr == "Unknown", "Unknown", "Tested"),
    tx_1L_grp  = case_when(
      str_detect(treatment_1L, "EGFR")          ~ "Targeted: EGFR TKI",
      str_detect(treatment_1L, "ALK")            ~ "Targeted: ALK Inhibitor",
      str_detect(treatment_1L, "Pembrolizumab")  ~ "Immunotherapy Mono",
      str_detect(treatment_1L, "Chemo-Immuno")   ~ "Chemo-Immunotherapy",
      TRUE                                       ~ "Chemotherapy Alone"
    )
  )

message(sprintf("Cohort: N = %d patients", nrow(df)))

# ==============================================================================
# TABLE 1 — Baseline Characteristics
# ==============================================================================
tab1_vars <- c("age","age_cat","sex","race","insurance","region",
               "stage","histology","ecog","smoking",
               "egfr_mut","alk_fusion","pdl1_expr",
               "cci","copd","diabetes","cvd")

tab1 <- CreateTableOne(
  vars      = tab1_vars,
  strata    = "tx_1L_grp",
  data      = df,
  factorVars = setdiff(tab1_vars, c("age","cci"))
)

tab1_mat <- print(tab1, showAllLevels = TRUE, quote = FALSE,
                  noSpaces = TRUE, printToggle = FALSE)

write.csv(tab1_mat,
          "study1_treatment_patterns/outputs/tables/table1_baseline.csv")
message("Table 1 saved.")

# ==============================================================================
# FIGURE 1 — First-Line Treatment Distribution Over Time
# ==============================================================================
tx_by_year <- df %>%
  group_by(year_dx, tx_1L_grp) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(year_dx) %>%
  mutate(pct = n / sum(n) * 100)

pal_tx <- c(
  "Targeted: EGFR TKI"      = "#2D6A9F",
  "Targeted: ALK Inhibitor"  = "#1B998B",
  "Immunotherapy Mono"       = "#E84855",
  "Chemo-Immunotherapy"      = "#FF9F1C",
  "Chemotherapy Alone"       = "#8B8C89"
)

fig1 <- ggplot(tx_by_year, aes(x = year_dx, y = pct, fill = tx_1L_grp)) +
  geom_area(alpha = 0.9, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = pal_tx, name = "1L Treatment") +
  scale_x_continuous(breaks = 2018:2023) +
  scale_y_continuous(labels = label_percent(scale = 1),
                     expand = c(0, 0)) +
  labs(
    title    = "Shift in First-Line NSCLC Treatment, 2018–2023",
    subtitle = "Increasing adoption of immunotherapy-based regimens over time",
    x        = "Year of Diagnosis",
    y        = "% of Patients",
    caption  = "Synthetic RWE data; portfolio demonstration only."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("study1_treatment_patterns/outputs/figures/fig1_tx_trends.png",
       fig1, width = 9, height = 5.5, dpi = 300)
message("Figure 1 saved.")

# ==============================================================================
# FIGURE 2 — Treatment Sequence Sankey / Alluvial
# ==============================================================================
seq_df <- df %>%
  filter(received_2L == 1, !is.na(treatment_2L)) %>%
  mutate(
    tx2_grp = case_when(
      treatment_2L == "Docetaxel"                ~ "Docetaxel",
      treatment_2L == "Docetaxel + Ramucirumab"  ~ "Docetaxel + Ramucirumab",
      treatment_2L == "Clinical Trial"           ~ "Clinical Trial",
      TRUE                                       ~ "Best Supportive Care"
    )
  ) %>%
  count(tx_1L_grp, tx2_grp) %>%
  rename(first_line = tx_1L_grp, second_line = tx2_grp, freq = n)

fig2 <- ggplot(seq_df,
               aes(axis1 = first_line, axis2 = second_line, y = freq)) +
  geom_alluvium(aes(fill = first_line), alpha = 0.75, width = 1/12) +
  geom_stratum(width = 1/5, fill = "grey92", colour = "grey50") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3.2, fontface = "bold") +
  scale_fill_manual(values = pal_tx, guide = "none") +
  scale_x_discrete(limits = c("1L Treatment","2L Treatment"),
                   expand = c(0.05, 0.05)) +
  labs(
    title    = "Treatment Sequencing: First-Line to Second-Line Therapy",
    subtitle = "Patient flow from 1L regimen to subsequent therapy",
    y        = "Number of Patients",
    caption  = "Restricted to patients who received 2L therapy."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.x = element_blank()
  )

ggsave("study1_treatment_patterns/outputs/figures/fig2_tx_sequence.png",
       fig2, width = 10, height = 6.5, dpi = 300)
message("Figure 2 saved.")

# ==============================================================================
# FIGURE 3 — Time to Treatment Initiation by Insurance Type
# ==============================================================================
tti_df <- df %>%
  mutate(tti_cat = cut(days_to_tx, c(0,14,30,60,Inf),
                       labels = c("≤14 days","15–30 days",
                                  "31–60 days",">60 days")))

tti_plot <- df %>%
  ggplot(aes(x = days_to_tx, fill = insurance)) +
  geom_histogram(bins = 40, alpha = 0.85, colour = "white", linewidth = 0.2) +
  facet_wrap(~ insurance, ncol = 2, scales = "free_y") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(
    title    = "Time from Diagnosis to First-Line Treatment Initiation",
    subtitle = "Stratified by insurance type — potential access disparities",
    x        = "Days to Treatment",
    y        = "Number of Patients",
    caption  = "Vertical dashed line = 30-day benchmark."
  ) +
  geom_vline(xintercept = 30, linetype = "dashed", colour = "grey30") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("study1_treatment_patterns/outputs/figures/fig3_time_to_tx.png",
       tti_plot, width = 9, height = 6, dpi = 300)
message("Figure 3 saved.")

# ==============================================================================
# ANALYSIS — Predictors of Biomarker Testing (Logistic Regression)
# ==============================================================================
df <- df %>%
  mutate(pdl1_tested = as.integer(pdl1_expr != "Unknown"),
         egfr_tested = as.integer(egfr_mut   != "Unknown"))

# PDL1 testing model
mod_pdl1 <- glm(pdl1_tested ~ age_cat + sex + race + insurance +
                  stage_grp + histology + ecog + year_dx,
                data   = df,
                family = binomial(link = "logit"))

or_pdl1 <- broom::tidy(mod_pdl1, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ round(.x, 3)))

write.csv(or_pdl1,
          "study1_treatment_patterns/outputs/tables/logistic_pdl1_testing.csv",
          row.names = FALSE)

# Forest plot of ORs
fig4 <- or_pdl1 %>%
  filter(!str_detect(term, "age_cat|year")) %>%
  mutate(term = str_replace_all(term, c(
    "sex"       = "Sex: ",
    "race"      = "Race: ",
    "insurance" = "Insurance: ",
    "stage_grp" = "Stage: ",
    "histology" = "Histology: ",
    "ecog"      = "ECOG: "
  ))) %>%
  ggplot(aes(x = estimate, y = fct_reorder(term, estimate))) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, colour = "#2D6A9F", linewidth = 0.8) +
  geom_point(size = 3, colour = "#E84855") +
  scale_x_log10() +
  labs(
    title    = "Predictors of PDL1 Biomarker Testing",
    subtitle = "Odds Ratios from multivariable logistic regression (log scale)",
    x        = "Odds Ratio (95% CI)",
    y        = NULL,
    caption  = "Reference categories: Female, White, Commercial insurance, Stage III, Adenocarcinoma."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("study1_treatment_patterns/outputs/figures/fig4_pdl1_forest.png",
       fig4, width = 9, height = 6, dpi = 300)
message("Figure 4 saved.")
message("\n✓ Study 1 analysis complete. All outputs written.")
