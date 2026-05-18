# =============================================================================
# Exercise 3: Dynamic Treatment Effects of Coethnicity on Road Expenditure
# Dataset: Burgess et al. (2015) — Kenya District-Year Panel 1963–2011
# Estimator: did_multiplegt_dyn (de Chaisemartin & D'Haultfoeuille 2024)
# =============================================================================
# NOTE: This script uses did_multiplegt_dyn to revisit the Burgess et al. (2015)
# setting. It does NOT replicate their main interacted TWFE specification
# (coethnic + coethnic × democracy). It estimates pooled dynamic effects of
# coethnicity and adds robustness and switchers-out exercises per the assignment.
# =============================================================================

# Set working directory to the ex_3 folder (adjust path if needed)
setwd("C:/Users/briac/Econ_Code/briac_s_ra_chaisemartin/ex_3")

# FIRST RUN: uncomment to install required packages
# install.packages(c("haven", "DIDmultiplegtDYN", "ggplot2", "dplyr"))
# DIDmultiplegtDYN 2.3.3 requires polars from r-universe (NOT on CRAN).
# Run this once when the network is available, then comment it out again:
# install.packages('polars', repos = 'https://rpolars.r-universe.dev')

library(haven)
library(polars)       # required by DIDmultiplegtDYN 2.3.3
library(DIDmultiplegtDYN)
library(ggplot2)
library(dplyr)

# Create output directories
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 1: DATA LOADING & PRE-ESTIMATION DIAGNOSTICS
# =============================================================================

df <- read_dta("exercise3_data.dta")

cat("\n=== DATA OVERVIEW ===\n")
str(df)
cat("\n=== SUMMARY OF KEY VARIABLES ===\n")
print(summary(df[, c("distnum", "year", "exp_dens_share", "president")]))
cat("\nAll variable names:", paste(names(df), collapse = ", "), "\n")

# --- 1.1 Panel Balance Check ---
# The package assumes evenly spaced time periods. A year missing for ALL groups
# must be filled before estimation.
cat("\n=== 1.1 PANEL BALANCE CHECK ===\n")
n_districts <- length(unique(df$distnum))
n_years     <- length(unique(df$year))
cat(sprintf("Districts: %d | Years: %d | Expected obs: %d | Actual obs: %d\n",
            n_districts, n_years, n_districts * n_years, nrow(df)))

year_counts <- df %>% group_by(year) %>% summarise(n_obs = n(), .groups = "drop")
missing_years <- year_counts %>% filter(n_obs < n_districts)
if (nrow(missing_years) > 0) {
  cat("WARNING: Years with fewer than expected districts:\n")
  print(as.data.frame(missing_years))
} else {
  cat("Panel is balanced: every year has", n_districts, "district observations.\n")
}

# --- 1.2 Treatment Switching Diagnostic ---
cat("\n=== 1.2 TREATMENT SWITCHING DIAGNOSTIC ===\n")
switch_diag <- df %>%
  arrange(distnum, year) %>%
  group_by(distnum) %>%
  mutate(
    switch_in  = president == 1 & lag(president, default = first(president)) == 0,
    switch_out = president == 0 & lag(president, default = first(president)) == 1
  ) %>%
  summarise(
    n_in          = sum(switch_in,  na.rm = TRUE),
    n_out         = sum(switch_out, na.rm = TRUE),
    always_zero   = all(president == 0),
    .groups = "drop"
  )

cat(sprintf("Districts with ≥1 switch IN  (0→1): %d\n", sum(switch_diag$n_in  > 0)))
cat(sprintf("Districts with ≥1 switch OUT (1→0): %d  [key for Q3]\n", sum(switch_diag$n_out > 0)))
cat(sprintf("Districts switching more than once (either direction): %d\n",
            sum((switch_diag$n_in + switch_diag$n_out) > 1)))
cat(sprintf("Never-treated districts (president=0 throughout): %d\n",
            sum(switch_diag$always_zero)))

# First-switch-in timing (which years did districts first become coethnic?)
first_switch_in_timing <- df %>%
  arrange(distnum, year) %>%
  group_by(distnum) %>%
  mutate(switch_in = president == 1 & lag(president, default = 0) == 0) %>%
  filter(switch_in) %>%
  slice(1) %>%   # first switch per district
  ungroup() %>%
  group_by(year) %>%
  summarise(n_districts = n(), .groups = "drop") %>%
  arrange(year)

cat("\nFirst switch-in timing (districts gaining coethnicity, by year):\n")
print(as.data.frame(first_switch_in_timing))

# --- 1.3 Available Economic-Geography Variables (for Q2 Column 4) ---
# Burgess et al. Table 1 Column 4 adds main-highway, border-district, and
# distance-to-Nairobi interactions. Check whether these are in the dataset.
cat("\n=== 1.3 ECONOMIC-GEOGRAPHY VARIABLE CHECK (for Q2 Col 4) ===\n")
geo_candidates <- c("highway", "main_road", "border", "dist_nairobi",
                    "dist_capital", "road_dist", "border_district")
found_geo <- names(df)[tolower(names(df)) %in% tolower(geo_candidates)]
if (length(found_geo) == 0) {
  cat("Economic-geography variables NOT found in dataset.\n")
  cat("Column 4 will match Column 3 (all available economic controls).\n")
  cat("This limitation is documented in the report.\n")
} else {
  cat("Found:", paste(found_geo, collapse = ", "), "\n")
}

# --- 1.4 Build Q2 Interaction Terms (baseline characteristic × time trend) ---
# These mimic θ(X_{d,1963} × [t - 1963]) from the benchmark regression.
# Time-varying by construction; the identifying assumption holds trivially.
df <- df %>% mutate(
  pop_trend       = pop1962         * (year - 1963),
  area_trend      = area            * (year - 1963),
  urbrate_trend   = urbrate1962     * (year - 1963),
  earnings_trend  = earnings        * (year - 1963),
  wage_trend      = wage_employment * (year - 1963),
  cashcrops_trend = value_cashcrops * (year - 1963)
)

trend_vars <- c("pop_trend", "area_trend", "urbrate_trend",
                "earnings_trend", "wage_trend", "cashcrops_trend")
na_counts <- sapply(trend_vars, function(v) sum(is.na(df[[v]])))
cat("\n=== 1.4 MISSING VALUES IN TREND VARIABLES ===\n")
print(na_counts)

# --- 1.5 Q2 Feasibility: Control Cell Counts by Baseline Treatment ---
# The package residualizes controls using only control (g,t) cells within each
# baseline-treatment group. If the number of controls is large relative to
# available cells, the package may warn about overfitting.
cat("\n=== 1.5 Q2 FEASIBILITY: CONTROL CELLS BY BASELINE TREATMENT ===\n")
baseline_treat <- df %>%
  arrange(distnum, year) %>%
  group_by(distnum) %>%
  summarise(baseline_pres = first(president), .groups = "drop")
df <- left_join(df, baseline_treat, by = "distnum")

cell_counts <- df %>%
  group_by(baseline_pres) %>%
  summarise(
    n_obs      = n(),
    n_districts = n_distinct(distnum),
    .groups = "drop"
  )
cat("Observations and districts by baseline treatment status:\n")
print(as.data.frame(cell_counts))
cat(sprintf("Col 3 uses 6 controls; monitor package warnings for overfitting.\n"))

# =============================================================================
# SECTION 2: Q1 — MAIN DYNAMIC EFFECTS (POOLED COETHNICITY EFFECT)
# =============================================================================
# ESTIMAND NOTE: did_multiplegt_dyn estimates horizon-specific ATTs by comparing
# switchers to groups with the same period-one treatment whose treatment has not
# yet changed at the relevant horizon. This is a POOLED effect of coethnicity,
# averaging over autocratic and democratic episodes. It differs from Burgess's
# main Table 1 estimand, which includes a coethnic × democracy interaction term
# separating effects by political regime. The pooled dynamic path may compress
# the ATT by mixing high-effect autocratic spells with lower-effect democratic
# ones. See report Section 3 for discussion.
#
# CHOICES:
#   effects = 7       : covers ~7 years post-switch; within the minimum
#                       presidential term in the sample (Kibaki: 9 years).
#   placebo = 4       : standard; sufficient to assess pre-trends (must ≤ effects).
#   cluster = distnum : treatment assigned at district level; specifying this
#                       makes the clustering choice explicit (package clusters
#                       at group level by default).
#   same_switchers = TRUE     : same composition of switchers used across all
#   same_switchers_pl = TRUE    effect and placebo horizons (prevents drift).
#   save_results      : writes estimates/SEs/CIs to CSV (print/summary have
#                       no return value; save_results is the reliable output path).
# =============================================================================

cat("\n\n=== Q1: MAIN DYNAMIC EFFECTS ===\n")

result_q1 <- did_multiplegt_dyn(
  df            = df,
  outcome       = "exp_dens_share",
  group         = "distnum",
  time          = "year",
  treatment     = "president",
  effects       = 7,
  placebo       = 4,
  cluster       = "distnum",
  same_switchers    = TRUE,
  same_switchers_pl = TRUE,
  # design option omitted: requires c(arg1, arg2) syntax in this package version,
  # not documented for binary-boolean use. Omitting does not affect estimates.
  save_results  = "output/tables/q1_raw.csv"
)

print(result_q1)

# plot.did_multiplegt_dyn renders a ggplot via print() then errors on a base-R call.
# tryCatch swallows that error; last_plot() captures what was just rendered.
tryCatch(plot(result_q1), error = function(e) invisible(NULL))
ggsave("output/event_study_q1.png", plot = ggplot2::last_plot(), width = 8, height = 5, dpi = 150, bg = "white")
cat("Saved: output/event_study_q1.png\n")

if (file.exists("output/tables/q1_raw.csv")) {
  q1_tab <- read.csv("output/tables/q1_raw.csv")
  cat("\nQ1 coefficient table (from save_results):\n")
  print(q1_tab)
} else {
  cat("Note: q1_raw.csv not created. Check that 'save_results' is a valid argument\n")
  cat("in your installed version of DIDmultiplegtDYN.\n")
}

# --- Q1 Sensitivity: Never-Switchers Only as Controls ---
# Restricts the control group to districts that never change treatment status.
# Provides a robustness check independent of not-yet-switched groups as controls.
cat("\n--- Q1 Sensitivity: only_never_switchers = TRUE ---\n")
n_never_sw <- sum(switch_diag$always_zero)
cat(sprintf("Never-switcher districts available: %d\n", n_never_sw))

if (n_never_sw >= 5) {
  result_q1_ns <- did_multiplegt_dyn(
    df        = df,
    outcome   = "exp_dens_share",
    group     = "distnum",
    time      = "year",
    treatment = "president",
    effects   = 7,
    placebo   = 4,
    cluster   = "distnum",
    only_never_switchers = TRUE,
    save_results = "output/tables/q1_ns_raw.csv"
  )
  print(result_q1_ns)
  tryCatch(plot(result_q1_ns), error = function(e) invisible(NULL))
  ggsave("output/event_study_q1_ns.png", plot = ggplot2::last_plot(), width = 8, height = 5, dpi = 150, bg = "white")
  cat("Saved: output/event_study_q1_ns.png\n")
} else {
  cat("Insufficient never-switchers for a powered sensitivity run. Skipping.\n")
}

# =============================================================================
# SECTION 3: Q2 — ROBUSTNESS: SEQUENTIAL CONTROLS (BURGESS TABLE 1, COLS 2–4)
# =============================================================================
# FEASIBILITY: Yes. did_multiplegt_dyn accepts time-varying controls via the
# 'controls' argument. The estimator first-differences controls and residualizes
# them against time FE, using only control (g,t) cells within each baseline-
# treatment group. Baseline-characteristic × linear-trend interactions are
# time-varying by construction, so this strategy is aligned with the estimator.
#
# CAVEAT: With 6 controls (Col 3), the package may warn about overfitting if
# the number of controls is large relative to control cells per baseline-
# treatment group. Record any such warnings.
#
# SCOPE LIMITATION: If economic-geography variables (main-highway, border,
# distance to Nairobi) are absent from the dataset, Column 4 repeats Column 3.
# This partial non-replication is documented in the report.
# =============================================================================

cat("\n\n=== Q2: SEQUENTIAL ROBUSTNESS CHECKS ===\n")

# Column 2: Population + Area + Urbanization × trend
cat("\n--- Q2 Column 2: Demographics × trend ---\n")
result_q2_col2 <- did_multiplegt_dyn(
  df        = df,
  outcome   = "exp_dens_share",
  group     = "distnum",
  time      = "year",
  treatment = "president",
  effects   = 7, placebo = 4, cluster = "distnum",
  # same_switchers omitted for Q2: combining it with controls on the small
  # baseline=1 group (7 districts) causes a polars Rust panic.
  controls  = c("pop_trend", "area_trend", "urbrate_trend"),
  save_results = "output/tables/q2_col2_raw.csv"
)
print(result_q2_col2)
tryCatch(plot(result_q2_col2), error = function(e) invisible(NULL))
ggsave("output/event_study_q2_col2.png", plot = ggplot2::last_plot(), width = 8, height = 5, dpi = 150, bg = "white")
cat("Saved: output/event_study_q2_col2.png\n")

# Column 3: + Earnings + Wage Employment + Cash Crops × trend
cat("\n--- Q2 Column 3: + Economic controls × trend ---\n")
result_q2_col3 <- did_multiplegt_dyn(
  df        = df,
  outcome   = "exp_dens_share",
  group     = "distnum",
  time      = "year",
  treatment = "president",
  effects   = 7, placebo = 4, cluster = "distnum",
  controls  = c("pop_trend", "area_trend", "urbrate_trend",
                "earnings_trend", "wage_trend", "cashcrops_trend"),
  save_results = "output/tables/q2_col3_raw.csv"
)
print(result_q2_col3)
tryCatch(plot(result_q2_col3), error = function(e) invisible(NULL))
ggsave("output/event_study_q2_col3.png", plot = ggplot2::last_plot(), width = 8, height = 5, dpi = 150, bg = "white")
cat("Saved: output/event_study_q2_col3.png\n")

# Column 4: + Economic-geography × trend
# Burgess Table 1 Col 4 adds main-highway, border-district, and distance-to-
# Nairobi interactions. If these variables are absent (see Section 1.3),
# uncomment the additional controls below and adjust variable names accordingly.
cat("\n--- Q2 Column 4: + Economic-geography × trend (if available) ---\n")
col4_controls <- c("pop_trend", "area_trend", "urbrate_trend",
                   "earnings_trend", "wage_trend", "cashcrops_trend")
# If economic-geography variables exist in the dataset, build trend interactions
# (e.g., highway_trend = highway * (year - 1963)) and add them here:
# col4_controls <- c(col4_controls, "highway_trend", "border_trend", "dist_nairobi_trend")

result_q2_col4 <- did_multiplegt_dyn(
  df        = df,
  outcome   = "exp_dens_share",
  group     = "distnum",
  time      = "year",
  treatment = "president",
  effects   = 7, placebo = 4, cluster = "distnum",
  controls  = col4_controls,
  save_results = "output/tables/q2_col4_raw.csv"
)
print(result_q2_col4)
tryCatch(plot(result_q2_col4), error = function(e) invisible(NULL))
ggsave("output/event_study_q2_col4.png", plot = ggplot2::last_plot(), width = 8, height = 5, dpi = 150, bg = "white")
cat("Saved: output/event_study_q2_col4.png\n")

# =============================================================================
# SECTION 4: Q3 — SWITCHERS OUT (EFFECTS OF LOSING COETHNICITY)
# =============================================================================
# Per Web Appendix 1.6 of de Chaisemartin & D'Haultfoeuille (2024):
# "Switchers out" are groups whose average treatment AFTER they switch is lower
# than their baseline treatment. In a binary treatment, this generally coincides
# with groups experiencing 1→0 transitions. However, multi-switcher districts
# (those experiencing both upward and downward treatment changes) may be dropped
# by default if the package detects ambiguous switching direction. The
# save_sample=TRUE argument tags each cell's role (switcher-in, switcher-out, or
# control) and records the effect number, allowing verification of the actual
# Q3 sample composition.
#
# SUBSTANTIVE INTERPRETATION: Positive Q1 coefficients indicate coethnicity
# increases road expenditure. Q3 tests whether this advantage is symmetrically
# reversed when a district loses coethnicity (new president from a different
# ethnic group), or whether built road infrastructure persists (hysteresis),
# in which case Q3 coefficients would be smaller in absolute value than Q1's
# or decay more slowly.
#
# HORIZON RULE: We report the maximum horizon at which the package identifies
# at least 10 switcher-out observations. This threshold is stated a priori.
# =============================================================================

cat("\n\n=== Q3: SWITCHERS OUT ===\n")
n_out_districts <- sum(switch_diag$n_out > 0)
cat(sprintf("Districts with ≥1 switch OUT (1→0): %d\n", n_out_districts))
if (n_out_districts < 5) {
  cat("WARNING: Very few switchers-out; Q3 may be underpowered.\n")
}

# Q3: SWITCHERS OUT — NOT FEASIBLE IN THIS DATASET
# -----------------------------------------------------------------------
# Attempting did_multiplegt_dyn with switchers = "out" raises:
#   "No treatment effect can be estimated. Design Restriction 1 in
#    de Chaisemartin & D'Haultfoeuille (2024) is not satisfied.
#    All groups experience their first treatment change at the same date."
#
# Root cause: The only potential switchers-out are the 7 Kenyatta/Kibaki
# Kikuyu districts (period-one treatment = 1). ALL seven lose coethnicity
# simultaneously in 1979 when Moi takes office. Design Restriction 1
# requires comparison units with the same period-one treatment (president=1)
# whose treatment has not yet changed at the relevant horizon. Because every
# baseline-treated district switches out in the same year, no valid
# comparison units exist at any horizon, and the estimator cannot be computed.
#
# The 6 Moi Kalenjin districts have period-one treatment = 0 and cannot
# serve as switchers-out (they are switchers-in relative to their baseline).
#
# Answer to Q3: Not feasible in this dataset. See report Section 4.3.

cat("\n\n=== Q3: SWITCHERS OUT ===\n")
cat("Diagnostic: checking switcher-out composition...\n")

# Count districts that have period-one treatment = 1 (potential switchers-out)
baseline1_districts <- sum(switch_diag$n_out > 0 &
  (df %>% arrange(distnum, year) %>% group_by(distnum) %>%
     summarise(b = first(president), .groups = "drop"))$b == 1)

# Check how many distinct years the period-one=1 group first switches out
first_switchout_timing <- df %>%
  arrange(distnum, year) %>%
  group_by(distnum) %>%
  mutate(switch_out = president == 0 & lag(president, default = first(president)) == 1,
         base_pres  = first(president)) %>%
  filter(switch_out, base_pres == 1) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(year) %>%
  summarise(n_districts = n(), .groups = "drop")

cat("Year(s) when baseline-treated districts first switch OUT:\n")
print(as.data.frame(first_switchout_timing))

cat("\nConclusion: All 7 baseline-coethnic (Kenyatta/Kibaki Kikuyu) districts\n")
cat("switch out simultaneously in 1979 when Moi takes office.\n")
cat("Design Restriction 1 of dCDH (2024) requires comparison units with the\n")
cat("same period-one treatment (president=1) whose treatment has not yet changed.\n")
cat("Since no such units exist, did_multiplegt_dyn with switchers='out' cannot\n")
cat("produce estimates. Q3 answer: NOT FEASIBLE in this dataset.\n")

# =============================================================================
# SECTION 5: SESSION INFO (for reproducibility)
# =============================================================================
cat("\n\n=== SESSION INFO ===\n")
sessionInfo()
