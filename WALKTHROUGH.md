# Portfolio Walkthrough Guide
## For Use in Interviews, Job Talks & Employer Presentations

This guide helps you walk any employer through your portfolio confidently — whether they are a pharmacoepidemiologist, a data scientist, or a hiring manager with a clinical background.

---

## How to Frame the Portfolio (Opening Statement)

> "I built this portfolio to demonstrate the three core skill sets I'd bring to a real-world evidence role: describing how drugs are actually used in practice, evaluating comparative effectiveness with proper confounding control, and identifying safety signals. I used both R and SAS because most pharma teams use one or both, and I wanted to show I'm comfortable in either environment. All analyses follow published RWE guidelines — the FDA RWE Framework, PCORnet methods, and RECORD-PE reporting standards."

---

## STUDY 1 WALKTHROUGH: Treatment Patterns in NSCLC

### What to say when asked "What did you do here?"

> "This study asks a common real-world question: how are patients actually being treated, and is it aligned with clinical guidelines? I built a cohort of 3,000 stage III/IV NSCLC patients and described their first-line treatment patterns from 2018 to 2023. The headline finding is the dramatic shift from chemotherapy to immunotherapy-based regimens — which mirrors what we've seen in published claims data studies."

### Walk through the figures

**Figure 1 — Stacked area chart (treatment over time)**
> "This is the core descriptive output. Each band represents a treatment class. You can see chemo-immunotherapy combinations dominating by 2021–2022, which tracks with the pivotal KEYNOTE trials entering clinical practice. EGFR TKIs and ALK inhibitors are flat narrow bands — as you'd expect, since only ~15–20% of patients are biomarker-positive."

**Figure 2 — Sankey/alluvial diagram (treatment sequencing)**
> "This is something I find really valuable in RWE — understanding the patient journey beyond a single treatment line. The alluvial plot shows where patients went after first-line therapy. You can see that Pembrolizumab patients commonly transition to docetaxel-based second-line, consistent with standard guidelines."

**Figure 3 — Time to treatment by insurance**
> "This one is more of a health disparities question. Medicaid patients had a longer median time to treatment initiation than commercially insured patients. In a real study, you'd want to adjust for clinical severity, distance to oncology centres, etc. But the raw pattern raises a legitimate equity question."

**Figure 4 — Forest plot: predictors of PDL1 testing**
> "This is the logistic regression output. I'm modelling who gets PDL1 biomarker testing, because without testing, patients can't access immunotherapy. The OR below 1 for squamous histology is interesting — it partly reflects historical prescribing before squamous-cell guidelines were updated."

### Technical depth — if asked to go deeper

- **Why logistic regression and not a log-binomial model?**  
  "PDL1 testing rates are around 70–80%, so the ORs from logistic regression won't closely approximate relative risks. A log-binomial model would give prevalence ratios directly, which are more interpretable in this context. In a real study I'd report both."

- **What variables would you add with real data?**  
  "Geographic distance to academic cancer centre, time period (important given guideline changes), and whether the patient was seen by a medical oncologist vs primary care first."

---

## STUDY 2 WALKTHROUGH: Survival Analysis in mCRC

### Opening

> "Study 2 is a comparative effectiveness question — does FOLFOX give better survival than FOLFIRI in metastatic colorectal cancer? This is perfect for demonstrating why you can't just run a Kaplan-Meier in RWE and call it a day. Treatment assignment isn't random: sicker patients with worse performance status and right-sided primaries were more likely to receive FOLFIRI. So I built a full methodological progression: crude → adjusted → IPTW weighted → landmark."

### Walking through the methods

**Kaplan-Meier**
> "The KM gives you the unadjusted picture. FOLFOX looks better, but we know the groups aren't balanced — so this is our starting point, not our answer."

**Cox PH — adjusted**
> "Multivariable Cox adds covariates as controls. The HR attenuates from 0.82 to ~0.86 once you adjust for ECOG, primary site, KRAS/BRAF status, and comorbidities. That attenuation is the confounding signal."

**IPTW — propensity score weighting**
> "IPTW is my preferred method for active-comparator studies. I estimate the propensity score — the probability of receiving FOLFOX given baseline covariates — and use inverse probability weights to create a pseudo-population where the treatment groups are balanced. The love plot shows the standardised mean differences before and after weighting. Convention is <0.10 is acceptable balance."

**Proportional hazards check**
> "This is something interviewers love to ask about. I test the PH assumption using Schoenfeld residuals — plotting them against time. If they're flat, the assumption holds. The formal test is a correlation between residuals and time. I also look for crossing survival curves in the KM."

**Landmark analysis**
> "The landmark analysis conditions on patients being alive at 12 months. This addresses a specific bias: if patients assigned to FOLFIRI are more likely to die early (before receiving 2L therapy), a simple Cox model conflates the treatment effect with that early mortality. The landmark removes that window."

**Subgroup forest plot**
> "This is for interaction/effect modification. I'm not doing formal interaction tests here, but you can see the HR is more favourable in KRAS wild-type patients, which is biologically plausible. In a real study, I'd include a p-for-interaction test."

### Technical questions

- **Why ATE and not ATT for IPTW?**  
  "ATE estimates the effect if the *entire* population received each treatment — a population-level policy question. ATT estimates the effect among those who actually received the treatment. For comparative effectiveness of two active drugs, ATE is usually more relevant. ATT would make sense if one arm is an active drug and one is no-treatment."

- **How do you handle extreme PS weights?**  
  "I trim at the 1st/99th percentiles (you can see this in the IPTW section of the SAS code) and report stabilised weights. I check the effective sample size after weighting. If weights are extreme, that often signals violations of positivity."

---

## STUDY 3 WALKTHROUGH: Drug Safety with CDK4/6 Inhibitors

### Opening

> "Study 3 is a pharmacovigilance and drug safety question — are CDK4/6 inhibitors associated with cardiovascular adverse events? Ribociclib in particular has a known QTc prolongation effect, so this is clinically motivated. I use three complementary approaches: incidence rate comparison, a new-user active-comparator cohort, and a negative control outcome."

### Walking through the design

**New-user active-comparator design**
> "This design choice matters a lot. 'New-user' means I restrict to patients initiating therapy — no prevalent users who've already tolerated the drug for months. This eliminates immortal time bias and depletion-of-susceptibles. 'Active-comparator' means I compare CDK4/6 inhibitors to endocrine monotherapy, not to 'no treatment' — which is a stronger design because the comparator group has a similar indication."

**Negative control outcome (UTI)**
> "This is one of my favourite tools in pharmacoepidemiology. A UTI has no plausible biological link to CDK4/6 inhibitors. So if I see a strong association there, it means there's systematic confounding or a data quality issue — not a real drug effect. The fact that my negative control shows IRR ~1.02 gives me confidence in the study validity."

**Disproportionality / ROR**
> "The reporting odds ratio mimics the signal detection methods used on spontaneous reporting databases like FAERS or VigiBase. It's not causal — it's signal detection. An ROR > 2 with tight CIs might trigger further investigation."

**AE heatmap**
> "The heatmap gives a multi-drug, multi-AE picture at once. You can immediately see that neutropenia is a class effect (all three CDK4/6 inhibitors show high rates) while the CV signal is more specific to ribociclib. That's a meaningful distinction for clinical decision-making."

### Technical questions

- **What's the difference between SCCS and a cohort design for safety?**  
  "Self-controlled case series uses only cases and uses each patient as their own control across time windows (pre-exposure, exposure, post-exposure). It eliminates all time-invariant confounding, which is powerful. But it requires the outcome not to affect exposure, and it's inappropriate for chronic outcomes or fatal events that censor follow-up. For a CV event that's relatively acute and doesn't eliminate follow-up immediately, it's a strong complement to the cohort design."

- **Why exact Poisson CIs for incidence rates?**  
  "Standard normal approximation CIs perform poorly with small event counts. Exact Poisson CIs (based on the chi-square distribution) are more reliable. In SAS this is the `CINV()` function."

---

## General Interview Tips

### Questions you might be asked and how to answer them

**"What's your favourite method in pharmacoepidemiology?"**
> "Propensity score methods — specifically IPTW — because they force you to explicitly model and communicate the source of confounding, rather than burying it in a covariate list. The love plot is a transparent way to show reviewers and regulators that you've achieved balance."

**"How would you handle unmeasured confounding?"**
> "Several ways. Negative control outcomes detect it empirically. Sensitivity analysis using E-values quantifies how strong unmeasured confounding would need to be to explain away the result. If the E-value is very large relative to what's plausible, the finding is more robust."

**"Why both R and SAS?"**
> "Many pharma companies — especially in regulatory-facing roles — still use SAS for validated submissions. CDISC standards and FDA-accepted outputs typically come from SAS. But R is increasingly dominant for exploratory analyses, visualisations, and methods development. Being fluent in both makes me useful in any team configuration."

**"What would you do differently with real data?"**
> "I'd validate diagnoses using positive predictive value studies of the ICD codes; I'd address the healthcare utilisation bias (people who never see a doctor can't be diagnosed or treated); I'd apply RECORD-PE or STaRT-RWE reporting checklists; and I'd engage clinical experts to review the plausibility of every covariate choice."

---

## Methodological Glossary (Quick Reference)

| Term | Plain-English Definition |
|------|--------------------------|
| **IPTW** | Re-weighting patients so treated/untreated groups look similar on confounders |
| **SMD** | Standardised Mean Difference — measures covariate imbalance; goal <0.10 |
| **SCCS** | Self-Controlled Case Series — each patient is their own control |
| **Landmark analysis** | Conditioning on survival to a time point to remove early-period bias |
| **Negative control** | Outcome with no biological link to drug — tests for hidden bias |
| **Active comparator** | Comparing two active therapies rather than drug vs. no treatment |
| **New-user design** | Restricting to patients newly initiating therapy (no prevalent users) |
| **E-value** | How strong unmeasured confounding must be to explain away an association |
| **PH assumption** | Cox model assumes the HR is constant over time — tested via Schoenfeld |
| **ROR** | Reporting Odds Ratio — disproportionality signal in spontaneous reports |

---

*Prepared as part of pharmacoepidemiology portfolio — [Your Name]*
