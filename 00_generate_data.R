# ==============================================================================
# Study 1: Real-World Treatment Patterns in Advanced NSCLC
# Script:   00_generate_data.R
# Purpose:  Generate synthetic patient-level RWD mimicking a claims/EHR dataset
# Author:   [Your Name]
# Note:     All data are entirely simulated for portfolio purposes only.
#           No real patient data are used.
# ==============================================================================

set.seed(2024)
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(tidyr)
})

N <- 3000

# ── 1. Demographics ────────────────────────────────────────────────────────────
patients <- tibble(
  patient_id  = sprintf("PT%05d", 1:N),
  age         = pmax(pmin(round(rnorm(N, 65, 10)), 90), 30),
  sex         = sample(c("Male","Female"), N, TRUE, c(0.55, 0.45)),
  race        = sample(c("White","Black","Hispanic","Asian","Other"),
                       N, TRUE, c(0.65, 0.13, 0.12, 0.07, 0.03)),
  insurance   = sample(c("Commercial","Medicare","Medicaid","Other"),
                       N, TRUE, c(0.35, 0.45, 0.15, 0.05)),
  region      = sample(c("Northeast","South","Midwest","West"),
                       N, TRUE, c(0.22, 0.38, 0.22, 0.18)),

  # ── 2. Clinical characteristics ──────────────────────────────────────────────
  stage       = sample(c("IIIA","IIIB","IVA","IVB"),
                       N, TRUE, c(0.12, 0.10, 0.42, 0.36)),
  histology   = sample(c("Adenocarcinoma","Squamous Cell","Large Cell","NOS"),
                       N, TRUE, c(0.50, 0.30, 0.10, 0.10)),
  ecog        = sample(0:3, N, TRUE, c(0.25, 0.40, 0.25, 0.10)),
  smoking     = sample(c("Current","Former","Never"),
                       N, TRUE, c(0.30, 0.50, 0.20)),

  # ── 3. Biomarkers (oncology-specific) ────────────────────────────────────────
  egfr_mut    = sample(c("Positive","Negative","Unknown"),
                       N, TRUE, c(0.15, 0.65, 0.20)),
  alk_fusion  = sample(c("Positive","Negative","Unknown"),
                       N, TRUE, c(0.05, 0.80, 0.15)),
  pdl1_expr   = sample(c("High (>=50%)","Low (1-49%)","Negative (<1%)","Unknown"),
                       N, TRUE, c(0.30, 0.25, 0.25, 0.20)),

  # ── 4. Comorbidities (Charlson-style) ────────────────────────────────────────
  cci         = pmax(rpois(N, 1.8), 0),
  copd        = sample(0:1, N, TRUE, c(0.70, 0.30)),
  diabetes    = sample(0:1, N, TRUE, c(0.75, 0.25)),
  cvd         = sample(0:1, N, TRUE, c(0.72, 0.28)),

  # ── 5. Index date (diagnosis) ────────────────────────────────────────────────
  index_date  = sample(seq(as.Date("2018-01-01"),
                           as.Date("2023-06-30"), by = "day"),
                       N, TRUE)
)

# ── 6. Assign first-line treatment (biomarker-driven, realistic) ──────────────
assign_1L <- function(egfr, alk, pdl1, hist) {
  if (egfr == "Positive")         return("EGFR TKI (Osimertinib)")
  if (alk  == "Positive")         return("ALK Inhibitor (Alectinib)")
  if (pdl1 == "High (>=50%)" & hist != "Squamous Cell")
                                  return("Pembrolizumab mono")
  if (hist == "Squamous Cell")    return("Chemo-Immunotherapy")
  return(sample(c("Chemo-Immunotherapy","Chemotherapy alone"),
                1, prob = c(0.65, 0.35)))
}

patients <- patients %>%
  rowwise() %>%
  mutate(
    treatment_1L = assign_1L(egfr_mut, alk_fusion, pdl1_expr, histology)
  ) %>%
  ungroup()

# ── 7. Time-to-event outcomes ─────────────────────────────────────────────────
# OS median varies by treatment (months, Weibull scale)
os_scale <- case_when(
  patients$treatment_1L == "EGFR TKI (Osimertinib)"  ~ 38,
  patients$treatment_1L == "ALK Inhibitor (Alectinib)"~ 42,
  patients$treatment_1L == "Pembrolizumab mono"        ~ 26,
  patients$treatment_1L == "Chemo-Immunotherapy"       ~ 22,
  TRUE                                                 ~ 14
)

patients <- patients %>%
  mutate(
    os_months     = round(pmin(rweibull(N, shape = 1.2, scale = os_scale), 60), 1),
    os_event      = as.integer(runif(N) > 0.20),     # 80% observed events
    pfs_months    = round(pmin(os_months * runif(N, 0.4, 0.85), os_months), 1),
    pfs_event     = as.integer(runif(N) > 0.10),
    ttnt_months   = round(pmax(pfs_months + runif(N, 0, 3), 0.5), 1),
    received_2L   = as.integer(os_months > ttnt_months),
    treatment_2L  = ifelse(
      received_2L == 1,
      sample(c("Docetaxel","Docetaxel + Ramucirumab","Clinical Trial",
               "Best Supportive Care"), N, TRUE, c(0.35, 0.30, 0.20, 0.15)),
      NA_character_
    ),
    # Time from diagnosis to first treatment (days)
    days_to_tx    = round(pmax(rnorm(N, 21, 12), 1)),
    start_date_1L = index_date + days_to_tx
  )

# ── 8. Save ───────────────────────────────────────────────────────────────────
write.csv(patients, "study1_treatment_patterns/data/nsclc_cohort.csv",
          row.names = FALSE)
message(sprintf("Dataset saved: %d patients, %d variables", nrow(patients), ncol(patients)))
message("Treatment distribution:")
print(table(patients$treatment_1L))
