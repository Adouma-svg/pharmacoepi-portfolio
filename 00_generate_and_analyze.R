# ==============================================================================
# Study 3: Drug Safety — Cardiovascular Risk with CDK4/6 Inhibitors
#           in Hormone Receptor+ Breast Cancer (Self-Controlled Case Series +
#           New-User Active-Comparator Design)
# Script:   00_generate_and_analyze.R
# Methods:  Incidence rate comparison, IRR (Poisson), SCCS, disproportionality
#           (ROR), new-user active-comparator cohort, negative control outcomes
# ==============================================================================

set.seed(2026)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(lubridate); library(purrr); library(scales)
  library(survival); library(broom); library(forcats)
  library(patchwork); library(RColorBrewer)
})

N <- 4000

# ==============================================================================
# PART 1 — GENERATE SYNTHETIC DATASET
# ==============================================================================

breast_ca <- tibble(
  patient_id   = sprintf("BC%05d", 1:N),
  age          = pmax(pmin(round(rnorm(N, 58, 12)), 85), 25),
  sex          = "Female",
  race         = sample(c("White","Black","Hispanic","Asian","Other"),
                        N, TRUE, c(0.67, 0.14, 0.11, 0.06, 0.02)),
  # Tumour characteristics
  er_status    = sample(c("Positive","Negative"), N, TRUE, c(0.92, 0.08)),
  her2_status  = "Negative",
  grade        = sample(1:3, N, TRUE, c(0.20, 0.45, 0.35)),
  nodal_status = sample(c("N0","N1","N2","N3"), N, TRUE, c(0.35,0.35,0.20,0.10)),
  stage        = sample(c("II","IIIA","IIIB","IV"), N, TRUE, c(0.25,0.25,0.15,0.35)),
  # Comorbidities
  cci          = pmax(rpois(N, 1.3), 0),
  htn          = sample(0:1, N, TRUE, c(0.55, 0.45)),
  diabetes     = sample(0:1, N, TRUE, c(0.78, 0.22)),
  prior_cvd    = sample(0:1, N, TRUE, c(0.85, 0.15)),
  # Drug assignment (CDK4/6i vs endocrine mono — confounded by indication)
  treatment = sample(
    c("Palbociclib","Ribociclib","Abemaciclib","Endocrine Monotherapy"),
    N, TRUE, c(0.30, 0.25, 0.20, 0.25)
  ),
  # Partner endocrine
  endocrine    = sample(c("Letrozole","Anastrozole","Fulvestrant"),
                        N, TRUE, c(0.50, 0.30, 0.20)),
  index_date   = sample(seq(as.Date("2016-01-01"),
                            as.Date("2023-06-01"), by="day"), N, TRUE)
)

# Drug class flag
breast_ca <- breast_ca %>%
  mutate(
    cdk46i      = as.integer(treatment != "Endocrine Monotherapy"),
    age_cat     = cut(age, c(24,49,59,69,85),
                      labels = c("<50","50-59","60-69","70+")),
    # Follow-up (varies by censoring mechanisms)
    fu_months   = round(pmin(rweibull(N, 1.4, 36), 60), 1),

    # Cardiovascular adverse events (primary outcome)
    # Ribociclib has known QTc/CV signal; model this
    cv_base_rate = 0.03 +
      0.02 * (treatment == "Ribociclib") +
      0.01 * (treatment == "Abemaciclib") +
      0.015 * (prior_cvd == 1) +
      0.008 * (htn == 1) +
      0.001 * pmax(age - 55, 0),
    cv_event    = as.integer(runif(N) < pmin(cv_base_rate * fu_months / 12, 0.90)),
    cv_time     = ifelse(cv_event == 1,
                         round(fu_months * runif(N, 0.1, 0.9), 1),
                         fu_months),

    # Negative control outcome (UTI — should NOT be associated with CDK4/6i)
    uti_event   = as.integer(runif(N) < 0.08),   # ~8% background
    uti_time    = ifelse(uti_event==1, round(fu_months*runif(N,0.1,0.9),1), fu_months),

    # Grade 3+ neutropenia (known CDK4/6i class effect)
    neutro_rate = 0.02 + 0.18 * cdk46i + 0.05*(treatment=="Palbociclib"),
    neutropenia = as.integer(runif(N) < pmin(neutro_rate, 0.99)),
    neutro_time = ifelse(neutropenia==1,
                         round(fu_months*runif(N,0.05,0.5),1), fu_months),

    # Person-years
    person_years = fu_months / 12
  )

write.csv(breast_ca, "study3_drug_safety/data/breast_cancer_cohort.csv",
          row.names=FALSE)
message(sprintf("Breast cancer safety dataset: N=%d", N))

# ==============================================================================
# PART 2 — INCIDENCE RATE ANALYSIS
# ==============================================================================

# ── 2A. CV Event Incidence Rates per 100 person-years ──────────────────────
ir_table <- breast_ca %>%
  group_by(treatment) %>%
  summarise(
    n         = n(),
    n_events  = sum(cv_event),
    py        = sum(person_years),
    ir        = n_events / py * 100,
    ir_lower  = qchisq(0.025, 2*n_events) / (2*py) * 100,
    ir_upper  = (qchisq(0.975, 2*(n_events+1))) / (2*py) * 100,
    .groups = "drop"
  ) %>%
  mutate(across(c(ir, ir_lower, ir_upper), ~round(.x, 2)))

write.csv(ir_table, "study3_drug_safety/outputs/tables/cv_incidence_rates.csv",
          row.names=FALSE)

fig_ir <- ir_table %>%
  mutate(treatment = fct_reorder(treatment, ir)) %>%
  ggplot(aes(x=ir, y=treatment, colour=treatment)) +
  geom_errorbarh(aes(xmin=ir_lower, xmax=ir_upper),
                 height=0.25, linewidth=1.2) +
  geom_point(size=5) +
  geom_vline(xintercept = ir_table$ir[ir_table$treatment=="Endocrine Monotherapy"],
             linetype="dashed", colour="grey40") +
  scale_colour_manual(
    values = c("Palbociclib"="#2D6A9F","Ribociclib"="#E84855",
               "Abemaciclib"="#FF9F1C","Endocrine Monotherapy"="#8B8C89"),
    guide = "none"
  ) +
  labs(
    title    = "Cardiovascular Event Incidence Rates by CDK4/6 Inhibitor",
    subtitle = "Rate per 100 person-years (95% exact CI) | Dashed = endocrine comparator",
    x        = "Incidence Rate per 100 Person-Years",
    y        = NULL
  ) +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold"))

ggsave("study3_drug_safety/outputs/figures/fig1_incidence_rates.png",
       fig_ir, width=9, height=5, dpi=300)

# ── 2B. Poisson Regression — Incidence Rate Ratios ───────────────────────────
poisson_fit <- glm(
  cv_event ~ treatment + age + htn + prior_cvd + cci +
    offset(log(pmax(person_years, 0.01))),
  data   = breast_ca,
  family = poisson(link="log")
)

irr_results <- tidy(poisson_fit, exponentiate=TRUE, conf.int=TRUE) %>%
  filter(str_detect(term,"treatment")) %>%
  mutate(across(c(estimate,conf.low,conf.high), ~round(.x,3)),
         comparator = "Endocrine Monotherapy (ref)")

write.csv(irr_results, "study3_drug_safety/outputs/tables/poisson_IRR.csv",
          row.names=FALSE)

# ── 2C. Disproportionality Analysis (Reporting Odds Ratio) ───────────────────
# Simulates spontaneous reporting signal analysis (like FAERS/VigiBase)
rpt_data <- breast_ca %>%
  mutate(
    a_cv    = cv_event * cdk46i,        # CV events in CDK4/6i
    b_cv    = cv_event * (1-cdk46i),    # CV events in control
    a_all   = cdk46i,                   # all CDK4/6i reports
    b_all   = 1-cdk46i                  # all comparator reports
  )

A <- sum(rpt_data$a_cv)  # drug + event
B <- sum(rpt_data$b_cv)  # comparator + event
C <- sum(rpt_data$a_all) - A  # drug + no event
D <- sum(rpt_data$b_all) - B  # comparator + no event

ror   <- (A / C) / (B / D)
ror_l <- exp(log(ror) - 1.96 * sqrt(1/A + 1/B + 1/C + 1/D))
ror_u <- exp(log(ror) + 1.96 * sqrt(1/A + 1/B + 1/C + 1/D))

message(sprintf("ROR (CDK4/6i vs control, CV events): %.3f (95%% CI: %.3f–%.3f)",
                ror, ror_l, ror_u))

# ── 2D. Adverse Event Profile — Heatmap ──────────────────────────────────────
ae_summary <- breast_ca %>%
  pivot_longer(cols=c(cv_event, neutropenia, uti_event),
               names_to="ae_type", values_to="event") %>%
  mutate(ae_label = recode(ae_type,
    cv_event   = "Cardiovascular Event",
    neutropenia = "Grade 3+ Neutropenia",
    uti_event  = "Urinary Tract Infection\n(Negative Control)"
  )) %>%
  group_by(treatment, ae_label) %>%
  summarise(pct = mean(event) * 100, .groups="drop")

fig_heatmap <- ae_summary %>%
  ggplot(aes(x=treatment, y=ae_label, fill=pct)) +
  geom_tile(colour="white", linewidth=0.8) +
  geom_text(aes(label=sprintf("%.1f%%", pct)),
            size=4, fontface="bold", colour="white") +
  scale_fill_gradientn(
    colours = c("#F7F9FC","#5B9BD5","#1B4F8A","#0D2137"),
    name    = "Event Rate (%)"
  ) +
  labs(
    title    = "Adverse Event Profile — CDK4/6 Inhibitors vs Endocrine Monotherapy",
    subtitle = "UTI is a negative control outcome (expected null association)",
    x        = NULL, y = NULL
  ) +
  theme_minimal(base_size=12) +
  theme(
    plot.title  = element_text(face="bold"),
    axis.text.x = element_text(angle=25, hjust=1)
  )

ggsave("study3_drug_safety/outputs/figures/fig2_ae_heatmap.png",
       fig_heatmap, width=9, height=5, dpi=300)

# ── 2E. New-User Active-Comparator Cohort — Survival Curves ──────────────────
ac_cohort <- breast_ca %>%
  filter(treatment %in% c("Ribociclib","Endocrine Monotherapy"))

km_cv <- survfit(Surv(cv_time, cv_event) ~ treatment, data=ac_cohort)

km_df <- broom::tidy(km_cv) %>%
  mutate(strata = str_remove(strata, "treatment="))

fig_km <- ggplot(km_df, aes(x=time, y=1-estimate, colour=strata, fill=strata)) +
  geom_step(linewidth=1.2) +
  geom_ribbon(aes(ymin=1-conf.high, ymax=1-conf.low), alpha=0.15, colour=NA) +
  scale_colour_manual(values=c("Ribociclib"="#E84855",
                                "Endocrine Monotherapy"="#8B8C89"),
                      name="Treatment") +
  scale_fill_manual(values=c("Ribociclib"="#E84855",
                              "Endocrine Monotherapy"="#8B8C89"),
                    guide="none") +
  scale_y_continuous(labels=label_percent(), limits=c(0,0.35)) +
  labs(
    title    = "Cumulative Incidence of CV Events: Ribociclib vs Endocrine Monotherapy",
    subtitle = "New-user active-comparator design | Shaded area = 95% CI",
    x        = "Follow-up (Months)",
    y        = "Cumulative Incidence",
    caption  = "Synthetic safety dataset; portfolio demonstration only."
  ) +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold"), legend.position="top")

ggsave("study3_drug_safety/outputs/figures/fig3_cumulative_incidence.png",
       fig_km, width=9, height=5.5, dpi=300)

message("\n✓ Study 3 analysis complete.")
