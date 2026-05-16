/*=============================================================================
  Study 2: Comparative Effectiveness — OS in Metastatic Colorectal Cancer
  Program:  01_survival_analysis.sas
  Methods:  Kaplan-Meier (LIFETEST), Cox PH (PHREG), Propensity Score IPTW,
            PH Assumption (Schoenfeld), Landmark Analysis
  Note:     SAS 9.4 / SAS Viya compatible.
=============================================================================*/

options nodate nonumber ls=200;
title;

libname study2 "./study2_survival_analysis/data";

/* ── 0. Import ───────────────────────────────────────────────────────────── */
proc import datafile="./study2_survival_analysis/data/mcrc_cohort.csv"
    out=study2.mcrc raw dbms=csv replace;
    getnames=yes;
run;

/* ── 1. Data preparation ─────────────────────────────────────────────────── */
data mcrc;
    set study2.mcrc_raw;
    /* Numeric flags */
    trt_folfox = (treatment = "FOLFOX");  /* reference = FOLFIRI */
    age_ge65   = (age >= 65);
    high_ecog  = (ecog >= 2);
    label
        os_months  = "Overall Survival (months)"
        pfs_months = "PFS (months)"
        os_event   = "OS Event (1=death)"
        trt_folfox = "Treatment: FOLFOX (1=yes)"
    ;
run;

/* ── 2. Kaplan-Meier — Overall Survival ──────────────────────────────────── */
title "Figure 1. Kaplan-Meier Estimates of Overall Survival";
title2 "FOLFOX vs FOLFIRI — Metastatic Colorectal Cancer";

ods graphics on / imagename="KM_OS" imagefmt=png width=700px height=500px;
proc lifetest data=mcrc
    plots=survival(atrisk cl test)
    timelist= 6 12 18 24 36 48
    outsurv=km_os_out;
    time      os_months * os_event(0);
    strata    treatment / test=logrank;
run;
ods graphics off;

/* Median OS by group */
title "Median OS by Treatment";
proc lifetest data=mcrc;
    time os_months * os_event(0);
    strata treatment;
run;
title;

/* ── 3. Cox PH — Unadjusted ──────────────────────────────────────────────── */
title "Table 2a. Unadjusted Cox PH — Overall Survival";
proc phreg data=mcrc;
    class treatment (ref="FOLFIRI");
    model os_months * os_event(0) = treatment / ties=efron rl;
run;

/* ── 4. Cox PH — Multivariable Adjusted ──────────────────────────────────── */
title "Table 2b. Multivariable-Adjusted Cox PH — Overall Survival";
proc phreg data=mcrc;
    class treatment   (ref="FOLFIRI")
          sex         (ref="Female")
          kras        (ref="Wild-type")
          braf        (ref="Wild-type")
          primary_site(ref="Left colon")
          / param=ref;
    model os_months * os_event(0) =
        treatment age sex ecog cci primary_site
        kras braf n_met_sites liver_met
        / ties=efron rl;
    hazardratio treatment / cl=pl;
    output out=cox_adj_out xbeta=linpred;
run;
title;

/* ── 5. Proportional Hazards Assumption (Schoenfeld Residuals) ───────────── */
title "PH Assumption Check — Schoenfeld Residuals";
proc phreg data=mcrc;
    class treatment (ref="FOLFIRI") / param=ref;
    model os_months * os_event(0) = treatment age ecog cci
                                    / ties=efron;
    output out=schoen_out schoenfeld=sch_treatment;
run;

/* Correlation of Schoenfeld residuals with time (p<0.05 = PH violated) */
proc corr data=schoen_out spearman;
    var sch_treatment;
    with os_months;
run;
title;

/* ── 6. Propensity Score Estimation ──────────────────────────────────────── */
title "Propensity Score Model — Treatment Assignment";
proc logistic data=mcrc descending;
    class sex kras braf primary_site / param=ref;
    model trt_folfox = age sex ecog cci primary_site
                       kras braf liver_met n_met_sites;
    output out=ps_out p=ps_hat;
run;

/* IPTW weights (ATE) */
data ps_weighted;
    set ps_out;
    if treatment="FOLFOX" then iptw = 1 / ps_hat;
    else                       iptw = 1 / (1 - ps_hat);
    /* Stabilized weights */
    p_trt = 0.50; /* approx marginal P(trt) */
    if treatment="FOLFOX" then sw_iptw = p_trt / ps_hat;
    else                       sw_iptw = (1-p_trt) / (1-ps_hat);
    /* Trim extreme weights at 1st/99th percentile */
run;

/* Check weight distribution */
title "IPTW Weight Distribution";
proc means data=ps_weighted n mean std min p1 p5 p95 p99 max;
    var iptw sw_iptw;
    class treatment;
run;

/* ── 7. Weighted Cox (IPTW) ──────────────────────────────────────────────── */
title "Table 3. IPTW-Weighted Cox PH — ATE Estimate";
proc phreg data=ps_weighted covs(aggregate);
    class treatment (ref="FOLFIRI") / param=ref;
    model os_months * os_event(0) = treatment / ties=efron rl;
    weight sw_iptw;
    id patient_id;
    hazardratio treatment / cl=pl;
run;
title;

/* ── 8. Landmark Analysis at 12 Months ───────────────────────────────────── */
title "Landmark Analysis — Patients Alive at 12 Months";
data landmark;
    set mcrc;
    where os_months >= 12;
    os_lm    = os_months - 12;
    event_lm = os_event;
run;

proc phreg data=landmark;
    class treatment (ref="FOLFIRI") / param=ref;
    model os_lm * event_lm(0) = treatment age ecog cci / ties=efron rl;
run;
title;

/* ── 9. Subgroup Analyses ────────────────────────────────────────────────── */
%macro cox_subgroup(label=, where_clause=);
    title "Subgroup: &label";
    proc phreg data=mcrc;
        where &where_clause;
        class treatment (ref="FOLFIRI") / param=ref;
        model os_months * os_event(0) = treatment / ties=efron rl;
    run;
%mend;

%cox_subgroup(label=Age lt 65,       where_clause=%str(age < 65));
%cox_subgroup(label=Age ge 65,       where_clause=%str(age >= 65));
%cox_subgroup(label=Male,            where_clause=%str(sex="Male"));
%cox_subgroup(label=Female,          where_clause=%str(sex="Female"));
%cox_subgroup(label=ECOG 0-1,        where_clause=%str(ecog <= 1));
%cox_subgroup(label=ECOG 2,          where_clause=%str(ecog = 2));
%cox_subgroup(label=KRAS Mutant,     where_clause=%str(kras="Mutant"));
%cox_subgroup(label=KRAS WT,         where_clause=%str(kras="Wild-type"));
%cox_subgroup(label=Liver Mets Yes,  where_clause=%str(liver_met=1));
%cox_subgroup(label=Liver Mets No,   where_clause=%str(liver_met=0));
title;

%put NOTE: Study 2 SAS survival analysis complete.;
