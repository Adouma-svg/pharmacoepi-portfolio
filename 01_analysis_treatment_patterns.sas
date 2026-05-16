/*=============================================================================
  Study 1: Real-World Treatment Patterns in Advanced NSCLC
  Program:  01_analysis_treatment_patterns.sas
  Purpose:  Replicate key analyses in SAS — Table 1, treatment frequencies,
            time to treatment, logistic regression for biomarker testing
  Author:   [Your Name]
  Note:     SAS 9.4 / SAS Viya compatible. Uses Base SAS + STAT procedures.
=============================================================================*/

options nodate nonumber ls=200 ps=max;
title;

/* ── 0. Libname & import ─────────────────────────────────────────────────── */
libname study1 "./study1_treatment_patterns/data";

proc import datafile="./study1_treatment_patterns/data/nsclc_cohort.csv"
    out     = study1.nsclc_raw
    dbms    = csv
    replace;
    getnames = yes;
run;

/* ── 1. Data preparation ─────────────────────────────────────────────────── */
data nsclc;
    set study1.nsclc_raw;

    /* Age categories */
    length age_cat $10;
    if      age <  55 then age_cat = "<55";
    else if age <= 64 then age_cat = "55-64";
    else if age <= 74 then age_cat = "65-74";
    else                   age_cat = "75+";

    /* Stage grouping */
    length stage_grp $12;
    if stage in ("IIIA","IIIB") then stage_grp = "Stage III";
    else                             stage_grp = "Stage IV";

    /* Year of diagnosis */
    year_dx = year(input(index_date, yymmdd10.));

    /* Simplified treatment group */
    length tx_1L_grp $35;
    if      index(treatment_1L,"EGFR")         > 0 then tx_1L_grp = "Targeted: EGFR TKI";
    else if index(treatment_1L,"ALK")          > 0 then tx_1L_grp = "Targeted: ALK Inhibitor";
    else if index(treatment_1L,"Pembrolizumab")> 0 then tx_1L_grp = "Immunotherapy Mono";
    else if index(treatment_1L,"Chemo-Immuno") > 0 then tx_1L_grp = "Chemo-Immunotherapy";
    else                                                tx_1L_grp = "Chemotherapy Alone";

    /* Biomarker testing flags */
    pdl1_tested = (pdl1_expr ne "Unknown");
    egfr_tested = (egfr_mut  ne "Unknown");

    /* Time-to-treatment category */
    length tti_cat $15;
    if      days_to_tx <=  14 then tti_cat = "<=14 days";
    else if days_to_tx <=  30 then tti_cat = "15-30 days";
    else if days_to_tx <=  60 then tti_cat = "31-60 days";
    else                           tti_cat = ">60 days";

    label
        age         = "Age at Diagnosis (years)"
        age_cat     = "Age Category"
        sex         = "Sex"
        race        = "Race/Ethnicity"
        insurance   = "Insurance Type"
        region      = "Geographic Region"
        stage       = "Disease Stage"
        histology   = "Histologic Subtype"
        ecog        = "ECOG Performance Status"
        smoking     = "Smoking Status"
        egfr_mut    = "EGFR Mutation Status"
        alk_fusion  = "ALK Fusion Status"
        pdl1_expr   = "PDL1 Expression"
        cci         = "Charlson Comorbidity Index"
        tx_1L_grp   = "First-Line Treatment Group"
        days_to_tx  = "Days from Diagnosis to Treatment"
        pdl1_tested = "PDL1 Testing Performed (0/1)"
    ;
run;

/* ── 2. Table 1 — Baseline Characteristics by Treatment Group ────────────── */
ods excel file="./study1_treatment_patterns/outputs/tables/table1_sas.xlsx"
    options(sheet_name="Table1" embedded_titles="yes");

title "Table 1. Baseline Characteristics by First-Line Treatment Group";
proc freq data=nsclc;
    tables (sex race insurance region stage histology smoking
            egfr_mut alk_fusion pdl1_expr) * tx_1L_grp
           / nocol nopercent chisq;
run;

title "Continuous Variables by Treatment Group";
proc means data=nsclc n mean std median q1 q3;
    class tx_1L_grp;
    var age cci days_to_tx;
run;

ods excel close;

/* ── 3. Treatment Frequency Table ────────────────────────────────────────── */
title "First-Line Treatment Distribution";
proc freq data=nsclc order=freq;
    tables tx_1L_grp / out=tx_freq;
run;

/* Trend over time */
title "First-Line Treatment by Year of Diagnosis";
proc freq data=nsclc;
    tables year_dx * tx_1L_grp / norow nocol;
run;
title;

/* ── 4. Time to Treatment Initiation ─────────────────────────────────────── */
title "Time to Treatment — Summary by Insurance Type";
proc means data=nsclc n mean std median q1 q3 min max;
    class insurance;
    var days_to_tx;
run;

title "Time to Treatment Category Distribution";
proc freq data=nsclc;
    tables insurance * tti_cat / norow nopercent;
run;
title;

/* ── 5. Logistic Regression — Predictors of PDL1 Testing ─────────────────── */
title "Table 3. Logistic Regression: Predictors of PDL1 Biomarker Testing";

/* Ensure reference categories */
proc logistic data=nsclc descending;
    class sex        (ref="Female")
          race       (ref="White")
          insurance  (ref="Commercial")
          age_cat    (ref="55-64")
          stage_grp  (ref="Stage III")
          histology  (ref="Adenocarcinoma")
          / param=ref;
    model pdl1_tested = age_cat sex race insurance stage_grp histology ecog year_dx
                       / clodds=pl risklimits;
    ods output OddsRatios = orout;
run;
title;

/* ── 6. Treatment Sequence (1L → 2L) ─────────────────────────────────────── */
data seq;
    set nsclc;
    where received_2L = 1 and treatment_2L ne "";
run;

title "Treatment Sequencing: Second-Line Therapy by First-Line Treatment";
proc freq data=seq;
    tables tx_1L_grp * treatment_2L / nocol nopercent;
run;
title;

/* ── 7. Export summary dataset for reporting ─────────────────────────────── */
proc export data=nsclc
    outfile="./study1_treatment_patterns/data/nsclc_analysis_ready.csv"
    dbms=csv replace;
run;

%put NOTE: Study 1 SAS analysis complete.;
