# ==============================================================================
# shared_utils/utils.R
# Reusable functions for pharmacoepidemiology analyses
# ==============================================================================

# ── Table 1 helpers ───────────────────────────────────────────────────────────

#' Format mean (SD) for continuous variables
fmt_mean_sd <- function(x, digits = 1) {
  sprintf("%.*f (%.*f)", digits, mean(x, na.rm=TRUE),
          digits, sd(x, na.rm=TRUE))
}

#' Format median [IQR]
fmt_median_iqr <- function(x, digits = 1) {
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm=TRUE)
  sprintf("%.*f [%.*f, %.*f]", digits, q[2], digits, q[1], digits, q[3])
}

#' Format n (%) for categorical variables
fmt_n_pct <- function(x, level = NULL) {
  if (!is.null(level)) x <- x == level
  n   <- sum(x, na.rm=TRUE)
  pct <- mean(x, na.rm=TRUE) * 100
  sprintf("%d (%.1f%%)", n, pct)
}

# ── Incidence rate with exact Poisson CI ─────────────────────────────────────

#' Calculate incidence rate and exact Poisson 95% CI
#' @param events integer — number of events
#' @param person_time numeric — total person-time
#' @param scale numeric — multiplier (default 100 for per 100 PY)
calc_ir <- function(events, person_time, scale = 100) {
  ir    <- events / person_time * scale
  lower <- qchisq(0.025, 2 * events) / (2 * person_time) * scale
  upper <- qchisq(0.975, 2 * (events + 1)) / (2 * person_time) * scale
  tibble::tibble(events, person_time, ir, ir_lower=lower, ir_upper=upper)
}

# ── PS / IPTW diagnostics ─────────────────────────────────────────────────────

#' Standardized Mean Difference (SMD) before/after weighting
smd <- function(x, trt, wt = NULL) {
  if (is.null(wt)) wt <- rep(1, length(x))
  x1 <- x[trt==1]; w1 <- wt[trt==1]
  x0 <- x[trt==0]; w0 <- wt[trt==0]
  m1 <- sum(x1*w1)/sum(w1); m0 <- sum(x0*w0)/sum(w0)
  v1 <- sum(w1*(x1-m1)^2)/sum(w1)
  v0 <- sum(w0*(x0-m0)^2)/sum(w0)
  (m1 - m0) / sqrt((v1 + v0) / 2)
}

# ── Survival helpers ──────────────────────────────────────────────────────────

#' Extract median survival with 95% CI from survfit object
get_median_os <- function(sf) {
  s <- summary(sf)$table
  if (is.matrix(s)) {
    data.frame(
      strata   = rownames(s),
      median   = s[,"median"],
      lower_95 = s[,"0.95LCL"],
      upper_95 = s[,"0.95UCL"]
    )
  } else {
    data.frame(median=s["median"], lower_95=s["0.95LCL"], upper_95=s["0.95UCL"])
  }
}

#' Restricted mean survival time (RMST) difference at tau
rmst_diff <- function(time, event, trt, tau) {
  fit1 <- survival::survfit(survival::Surv(time[trt==1], event[trt==1]) ~ 1)
  fit0 <- survival::survfit(survival::Surv(time[trt==0], event[trt==0]) ~ 1)
  rmst <- function(fit) {
    with(summary(fit), {
      t <- c(0, time[time <= tau]); s <- c(1, surv[time <= tau])
      sum(diff(t) * s[-length(s)])
    })
  }
  list(rmst_trt=rmst(fit1), rmst_ctrl=rmst(fit0),
       diff=rmst(fit1)-rmst(fit0))
}

# ── Output helpers ────────────────────────────────────────────────────────────

#' Save ggplot with consistent defaults
save_fig <- function(p, path, w=9, h=6, dpi=300) {
  ggplot2::ggsave(path, p, width=w, height=h, dpi=dpi, bg="white")
  message(sprintf("Saved: %s", path))
}

#' Round all numeric columns in a data frame
round_df <- function(df, digits=3) {
  df %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, digits)))
}
