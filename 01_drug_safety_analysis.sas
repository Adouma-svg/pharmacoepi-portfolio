/*=============================================================================
  Study 3: Drug Safety — Cardiovascular Risk with CDK4/6 Inhibitors
  Program:  01_drug_safety_analysis.sas
  Methods:  Incidence rates (exact CI), Poisson regression, IRR,
            New-user active-comparator Cox, Negative control outcome
=============================================================================*/

options nodate nonumber ls=200;
title;
libname study3 "./study3_drug_safety/data";

proc import datafile="./study3_drug_safety/data/breast_cancer_cohort.csv"
    out=study3.bc_raw dbms=csv replace;
    getnames=yes;
run;

/* ── 1. Data preparation ─────────────────────────────────────────────────── */
data bc;
    set study3.bc_raw;
    cdk46i      = (treatment ne "Endocrine Monotherapy");
    age_ge65    = (age >= 65);
    length age_cat $8;
    if      age <  50 then age_cat = "<50";
    else if age <= 59 then age_cat = "50-59";
    else if age <= 69 then age_cat = "60-69";
    else                   age_cat = "70+";
    label
        cv_event     = "CV Event (1=yes)"
        neutropenia  = "Gr3+ Neutropenia (1=yes)"
        uti_event    = "UTI - Negative Control (1=yes)"
        person_years = "Person-years at risk"
        cdk46i       = "CDK4/6 Inhibitor Class (1=yes)"
    ;
run;

/* ── 2. Descriptive: Adverse Event Rates ─────────────────────────────────── */
title "Table 1. Adverse Event Rates by Treatment Group";
proc means data=bc n sum mean std;
    class treatment;
    var cv_event neutropenia uti_event person_years;
run;
title;

/* ── 3. Incidence Rates with Exact Poisson CI ────────────────────────────── */
title "Cardiovascular Event Incidence Rates (per 100 PY)";
proc means data=bc n sum;
    class treatment;
    var cv_event person_years;
    output out=rate_dat sum(cv_event)=n_events sum(person_years)=total_py;
run;

data rate_final;
    set rate_dat;
    where _type_=1;
    ir        = n_events / total_py * 100;
    /* Exact Poisson 95% CI */
    ir_lower  = cinv(0.025, 2*n_events)     / (2*total_py) * 100;
    ir_upper  = cinv(0.975, 2*(n_events+1)) / (2*total_py) * 100;
    format ir ir_lower ir_upper 8.3;
run;

proc print data=rate_final label noobs;
    var treatment n_events total_py ir ir_lower ir_upper;
run;
title;

/* ── 4. Poisson Regression — Incidence Rate Ratios ──────────────────────── */
title "Table 2. Poisson Regression: CV Event Incidence Rate Ratios";
proc genmod data=bc;
    class treatment (ref="Endocrine Monotherapy") / param=ref;
    model cv_event = treatment age htn prior_cvd cci
        / dist=poisson link=log offset=log_py;
    /* Create offset */
    /* Note: add log_py to dataset first */
run;

/* Add log offset */
data bc2;
    set bc;
    log_py = log(max(person_years, 0.001));
run;

proc genmod data=bc2;
    class treatment (ref="Endocrine Monotherapy")
          age_cat   (ref="50-59")
          / param=ref;
    model cv_event = treatment age_cat htn prior_cvd cci
        / dist=poisson link=log offset=log_py;
    estimate "Palbociclib vs Control"       treatment 1 0 0 -1 / exp;
    estimate "Ribociclib vs Control"        treatment 0 1 0 -1 / exp;
    estimate "Abemaciclib vs Control"       treatment 0 0 1 -1 / exp;
    lsmeans treatment / ilink diff exp;
run;
title;

/* ── 5. Negative Control Outcome Analysis ────────────────────────────────── */
title "Negative Control: UTI Incidence Rate Ratios (CDK4/6i vs Endocrine)";
proc genmod data=bc2;
    class cdk46i / param=ref;
    model uti_event = cdk46i age cci
        / dist=poisson link=log offset=log_py;
    estimate "CDK4/6i vs Endocrine" cdk46i 1 / exp;
run;
/* Expected: IRR ~1.0 — validates no systematic bias for unrelated outcome */
title;

/* ── 6. New-User Active-Comparator Cohort ────────────────────────────────── */
data ac_cohort;
    set bc;
    where treatment in ("Ribociclib","Endocrine Monotherapy");
run;

title "New-User Active-Comparator: Ribociclib vs Endocrine Monotherapy";
title2 "Unadjusted Cox PH — CV Events";
proc phreg data=ac_cohort;
    class treatment (ref="Endocrine Monotherapy") / param=ref;
    model cv_time * cv_event(0) = treatment / ties=efron rl;
run;

title2 "Adjusted Cox PH — CV Events";
proc phreg data=ac_cohort;
    class treatment (ref="Endocrine Monotherapy")
          age_cat   (ref="50-59")
          / param=ref;
    model cv_time * cv_event(0) =
        treatment age_cat htn prior_cvd cci diabetes
        / ties=efron rl;
    hazardratio treatment / cl=pl;
run;
title;

/* ── 7. Kaplan-Meier Curves — CV Event Free Survival ─────────────────────── */
ods graphics on / imagename="KM_CV" imagefmt=png width=700px height=500px;
title "Figure 3. Cumulative Incidence of CV Events";
proc lifetest data=ac_cohort
    plots=cif method=km timelist=6 12 18 24 36;
    time cv_time * cv_event(0);
    strata treatment / test=logrank;
run;
ods graphics off;
title;

/* ── 8. Neutropenia — Class Effect Confirmation ──────────────────────────── */
title "CDK4/6 Class Effect: Grade 3+ Neutropenia by Agent";
proc freq data=bc;
    tables treatment * neutropenia / nocol nopercent riskdiff;
run;

proc logistic data=bc descending;
    class treatment (ref="Endocrine Monotherapy") / param=ref;
    model neutropenia = treatment age cci;
    oddsratio treatment;
run;
title;

%put NOTE: Study 3 SAS safety analysis complete.;
