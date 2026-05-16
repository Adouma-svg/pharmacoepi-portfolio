# Pharmacoepidemiology & Real-World Evidence Portfolio

**[Your Name]** | PhD Candidate in Pharmacoepidemiology  
📧 [your.email@university.edu] | 🌐 [portfolio-url] | 💼 [LinkedIn] | 🐙 [GitHub]

---

## Overview

This portfolio demonstrates applied pharmacoepidemiology and real-world evidence (RWE) methods using **R** and **SAS** — the two dominant platforms in pharma/biotech RWE teams. All datasets are fully synthetic (no patient data) and designed to mirror the structure of real-world clinical databases (claims, EHR, registry).

Three oncology case studies cover the core methodological pillars expected in RWE/data analytics roles:

| Study | Disease | Design | Key Methods |
|-------|---------|--------|-------------|
| [1](#study-1) | Advanced NSCLC | Descriptive RWE | Table 1, treatment trends, logistic regression, forest plot |
| [2](#study-2) | Metastatic CRC | Comparative Effectiveness | KM, Cox PH, IPTW, landmark analysis, subgroup forest |
| [3](#study-3) | HR+ Breast Cancer | Drug Safety | Incidence rates, Poisson/IRR, new-user active-comparator, negative control |

---

## Repository Structure

```
pharmacoepi_portfolio/
├── study1_treatment_patterns/       # NSCLC treatment patterns & disparities
│   ├── R/
│   │   ├── 00_generate_data.R       # Synthetic data generation
│   │   └── 01_analysis_treatment_patterns.R
│   ├── SAS/
│   │   └── 01_analysis_treatment_patterns.sas
│   ├── data/                        # Generated CSV datasets
│   └── outputs/
│       ├── figures/                 # PNG figures
│       └── tables/                  # CSV result tables
│
├── study2_survival_analysis/        # mCRC comparative effectiveness
│   ├── R/
│   │   └── 00_generate_and_analyze.R
│   ├── SAS/
│   │   └── 01_survival_analysis.sas
│   ├── data/
│   └── outputs/
│
├── study3_drug_safety/              # CDK4/6i cardiovascular safety
│   ├── R/
│   │   └── 00_generate_and_analyze.R
│   ├── SAS/
│   │   └── 01_drug_safety_analysis.sas
│   ├── data/
│   └── outputs/
│
├── shared_utils/
│   ├── utils.R                      # Reusable R functions
│   └── macros.sas                   # Reusable SAS macros
│
└── docs/
    ├── WALKTHROUGH.md               # Interview walkthrough guide
    └── METHODS_GLOSSARY.md          # Epidemiology methods reference
```

---

## Study 1: Real-World Treatment Patterns in Advanced NSCLC {#study-1}

### Research Question
How have first-line treatment patterns shifted in advanced non-small cell lung cancer (NSCLC) following the introduction of immunotherapy and targeted therapies? Are there disparities in biomarker testing rates?

### Dataset
- **N = 3,000** patients with stage III–IV NSCLC (2018–2023)
- Variables: demographics, insurance, biomarker status (EGFR/ALK/PDL1), comorbidities, treatment assignment, time-to-treatment, treatment sequencing

### Methods
- **Table 1** — Standardised baseline characteristics by treatment group (`tableone`, `gtsummary`)
- **Stacked area chart** — Temporal shift in 1L treatment mix (2018–2023)
- **Alluvial/Sankey diagram** — 1L → 2L treatment sequencing (`ggalluvial`)
- **Time-to-treatment histogram** — Stratified by insurance type (access disparities)
- **Logistic regression + forest plot** — Predictors of PDL1 biomarker testing

### Key Findings (Simulated)
- Immunotherapy-based regimens increased from ~30% (2018) to ~65% (2023) of 1L treatments
- Medicaid patients had longer time-to-treatment than commercially insured patients (median 28 vs 19 days)
- Older age and squamous histology were associated with lower PDL1 testing rates

### Files
- `R/00_generate_data.R` — Reproducible synthetic data
- `R/01_analysis_treatment_patterns.R` — Full analysis pipeline
- `SAS/01_analysis_treatment_patterns.sas` — SAS equivalent (PROC FREQ, PROC LOGISTIC)

---

## Study 2: Comparative Effectiveness — FOLFOX vs FOLFIRI in Metastatic CRC {#study-2}

### Research Question
Does FOLFOX provide superior overall survival compared to FOLFIRI as first-line chemotherapy in metastatic colorectal cancer (mCRC), after adjusting for confounding by indication?

### Dataset
- **N = 2,500** patients with mCRC receiving first-line chemotherapy
- Variables: demographics, ECOG status, molecular markers (KRAS/BRAF/MSI), metastasis pattern, treatment, OS/PFS outcomes

### Methods
| Method | Purpose |
|--------|---------|
| Kaplan-Meier | Unadjusted OS/PFS curves with risk tables |
| Log-rank test | Unadjusted OS comparison |
| Multivariable Cox PH | Adjusted HR controlling for confounders |
| **IPTW (ATE)** | Propensity score weighting to balance treatment groups |
| Love plot / SMD | Covariate balance assessment before/after IPTW |
| Schoenfeld residuals | Proportional hazards assumption check |
| **Landmark analysis** | Avoid immortal time bias — conditioned on 12-month survival |
| Subgroup forest plot | Treatment effect heterogeneity across key subgroups |

### Why IPTW Matters Here
Treatment assignment is *not* random in RWE — sicker patients (higher ECOG, right-sided tumours) were more likely to receive FOLFIRI. Naive comparisons are biased. IPTW reweights the population to mimic randomisation, providing a more valid estimate of the average treatment effect.

### Key Findings (Simulated)
- Crude HR: FOLFOX vs FOLFIRI = 0.82 (95% CI: 0.73–0.92)
- IPTW-adjusted HR: 0.86 (95% CI: 0.77–0.96) — attenuated after confounding control
- No PH assumption violation detected; landmark analysis consistent with primary result
- Effect modification by KRAS status (stronger benefit in wild-type patients)

### Files
- `R/00_generate_and_analyze.R` — Data generation + full survival analysis
- `SAS/01_survival_analysis.sas` — PROC LIFETEST, PROC PHREG, PS weighting

---

## Study 3: CDK4/6 Inhibitor Cardiovascular Safety in HR+ Breast Cancer {#study-3}

### Research Question
Are CDK4/6 inhibitors (palbociclib, ribociclib, abemaciclib) associated with an elevated risk of cardiovascular events compared to endocrine monotherapy in patients with hormone receptor-positive breast cancer?

### Dataset
- **N = 4,000** female patients with HR+/HER2− breast cancer
- Variables: tumour stage, molecular markers, comorbidities (CVD, HTN, diabetes), CDK4/6 inhibitor identity, CV events, Grade 3+ neutropenia, UTI (negative control)

### Methods
| Method | Purpose |
|--------|---------|
| Incidence rates (exact Poisson CI) | Baseline event rate quantification |
| Poisson regression (IRR) | Rate ratio adjusting for confounders |
| Reporting Odds Ratio (ROR) | Disproportionality signal detection (à la FAERS) |
| **New-user active-comparator design** | Minimise prevalent user & indication bias |
| **Negative control outcome (UTI)** | Detect unmeasured confounding / bias |
| Adverse event profile heatmap | Multi-AE multi-drug visual comparison |
| Cumulative incidence curves | Time-to-event visualisation |

### Why Negative Controls Matter
If UTI (a biologically unrelated outcome) showed a *positive* association with CDK4/6 inhibitors, it would suggest systematic confounding or data artefact. A null result for the negative control increases confidence that the CV signal is real.

### Key Findings (Simulated)
- Ribociclib showed the highest CV event rate (consistent with known QTc effects)
- Palbociclib showed Grade 3+ neutropenia rates ~3× higher than endocrine monotherapy
- UTI (negative control): IRR ~1.02 — no signal, supporting validity of study design
- In the new-user active-comparator cohort, ribociclib HR for CV events: 1.48 (1.12–1.96)

### Files
- `R/00_generate_and_analyze.R` — Data + safety analysis pipeline
- `SAS/01_drug_safety_analysis.sas` — PROC GENMOD (Poisson), PROC PHREG, PROC LIFETEST

---

## Technical Stack

| Tool | Version | Primary Use |
|------|---------|-------------|
| R | ≥ 4.3 | All analyses + visualisations |
| SAS | 9.4 / Viya | Parallel implementation |
| `survival` | CRAN | Cox PH, KM, Schoenfeld |
| `WeightIt` | CRAN | Propensity score IPTW |
| `survminer` | CRAN | Publication-quality KM plots |
| `gtsummary` | CRAN | Table 1 generation |
| `ggalluvial` | CRAN | Treatment sequence diagrams |
| `tableone` | CRAN | SMD-based balance tables |

### R Package Installation
```r
install.packages(c(
  "tidyverse", "survival", "survminer", "WeightIt", "cobalt",
  "tableone", "gtsummary", "ggalluvial", "broom", "patchwork",
  "scales", "lubridate", "forcats", "RColorBrewer", "purrr"
))
```

---

## Running the Analyses

```bash
# Clone the repository
git clone https://github.com/[yourusername]/pharmacoepi-portfolio.git
cd pharmacoepi-portfolio

# Run all R analyses in sequence
Rscript study1_treatment_patterns/R/00_generate_data.R
Rscript study1_treatment_patterns/R/01_analysis_treatment_patterns.R
Rscript study2_survival_analysis/R/00_generate_and_analyze.R
Rscript study3_drug_safety/R/00_generate_and_analyze.R
```

For SAS: open each `.sas` file in SAS EG or SAS Studio and submit from the portfolio root directory.

---

## Methodological Competencies Demonstrated

- ✅ Real-world evidence study design (descriptive, CER, safety)
- ✅ Claims/EHR data structure and variable construction
- ✅ Confounding control: multivariable regression, propensity score IPTW
- ✅ Time-to-event analysis: KM, Cox PH, competing risks concepts
- ✅ Bias assessment: PH assumption, landmark analysis, negative controls
- ✅ Health disparities analysis (insurance, race, access)
- ✅ Dual proficiency: R (tidyverse + survival) and SAS (Base + STAT)
- ✅ Reproducible, well-documented, version-controlled code

---

*All data in this portfolio are entirely synthetic, generated solely for demonstration purposes. No real patient information is used or implied.*
