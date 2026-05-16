/*=============================================================================
  shared_utils/macros.sas
  Reusable SAS macros for pharmacoepidemiology analyses
=============================================================================*/

/* ── 1. Table 1 — n(%) for categorical variable ──────────────────────────── */
%macro tab1_cat(data=, var=, group=, out=);
    proc freq data=&data noprint;
        tables &var * &group / outpct out=&out;
    run;
%mend tab1_cat;

/* ── 2. Incidence rate with exact Poisson CI ─────────────────────────────── */
%macro calc_ir(data=, event=, py=, group=, scale=100, out=ir_out);
    proc means data=&data noprint;
        class &group;
        var &event &py;
        output out=&out(where=(_type_=1))
            sum(&event)=n_events
            sum(&py)=total_py;
    run;

    data &out;
        set &out;
        ir       = n_events / total_py * &scale;
        ir_lower = cinv(0.025, 2*n_events)     / (2*total_py) * &scale;
        ir_upper = cinv(0.975, 2*(n_events+1)) / (2*total_py) * &scale;
        format ir ir_lower ir_upper 8.3;
        label
            ir       = "Incidence Rate per &scale PY"
            ir_lower = "Lower 95% CI"
            ir_upper = "Upper 95% CI";
    run;

    proc print data=&out noobs label;
        var &group n_events total_py ir ir_lower ir_upper;
    run;
%mend calc_ir;

/* Usage example:
   %calc_ir(data=bc, event=cv_event, py=person_years, group=treatment);
*/

/* ── 3. Cox subgroup loop ─────────────────────────────────────────────────── */
%macro cox_loop(data=, time=, event=, trt=, ref=, subvar=, subval=);
    /* Runs Cox PH in a data subset */
    data _sub_;
        set &data;
        where &subvar = "&subval";
    run;
    title "Subgroup: &subvar = &subval";
    proc phreg data=_sub_;
        class &trt (ref="&ref") / param=ref;
        model &time * &event(0) = &trt / ties=efron rl;
    run;
    proc datasets nolist; delete _sub_; run;
%mend cox_loop;

/* ── 4. Standardized Mean Difference ─────────────────────────────────────── */
%macro smd(data=, var=, trt=, out=smd_out);
    /* Computes SMD for a continuous/binary variable between trt groups */
    proc means data=&data noprint;
        class &trt;
        var &var;
        output out=_means_ mean=mu var=v n=n;
    run;

    data &out;
        set _means_(where=(&trt^=.));
        proc sort; by &trt; run;

    data &out;
        merge _means_(where=(&trt=1) rename=(mu=mu1 v=v1))
              _means_(where=(&trt=0) rename=(mu=mu0 v=v0));
        smd = (mu1 - mu0) / sqrt((v1 + v0) / 2);
        abs_smd = abs(smd);
        variable = "&var";
        keep variable smd abs_smd;
    run;

    proc datasets nolist; delete _means_; run;
%mend smd;

/* ── 5. Kaplan-Meier median + CI table ───────────────────────────────────── */
%macro km_table(data=, time=, event=, strata=);
    title "Kaplan-Meier Median Survival — &strata";
    proc lifetest data=&data;
        time &time * &event(0);
        strata &strata / test=logrank;
    run;
    title;
%mend km_table;

/* ── 6. Propensity score trimming (symmetric) ────────────────────────────── */
%macro trim_ps(data=, ps=ps_hat, pct_trim=0.01, out=trimmed);
    proc univariate data=&data noprint;
        var &ps;
        output out=_pcts_
            pctlpts  = %sysevalf(&pct_trim*100) %sysevalf(100-&pct_trim*100)
            pctlpre  = p;
    run;

    data _null_;
        set _pcts_;
        call symputx('ps_low',  p&pct_trim*100);
        call symputx('ps_high', p%sysevalf(100-&pct_trim*100));
    run;

    data &out;
        set &data;
        where &ps >= &ps_low and &ps <= &ps_high;
    run;
    %put NOTE: PS trimmed at [&ps_low, &ps_high]. N=;
    proc sql; select count(*) as n_after_trim from &out; quit;
%mend trim_ps;

%put NOTE: Pharmacoepi SAS macros loaded successfully.;
