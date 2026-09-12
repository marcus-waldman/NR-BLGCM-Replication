#' Step 4: build manuscript Figures 1, 2, 4 and 5 into figures/.
#'
#'   Fig 1  true-score vs observed trajectory (conceptual; no fit needed)
#'   Fig 2  prior / likelihood / posterior on residual SD (conceptual; no fit needed)
#'   Fig 4  36-week OGTT trajectories by arm + AUC treatment-effect insets
#'   Fig 5  prior-sensitivity forest plot of the adjusted AUC treatment effect
#'
#' Figure 3 (the path diagram) was drawn outside R and is not produced here.
#' Figure 6 is built by scripts/05_make_figure6_priorsense.R.
#'
#' Run from the repository root:  Rscript scripts/04_make_figures.R
#' Figures 4 and 5 need derived/estimands.rds (step 3).

local({
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a) == 1) {
    setwd(normalizePath(file.path(dirname(sub("^--file=", "", a)), "..")))
  }
})
suppressMessages({
  library(tidyverse)
  library(patchwork)
})
source("R/_setup.R")

#-------------------------------------------------------------------------------
# Figures 1 and 2: conceptual, seeded, independent of the data and fits
#-------------------------------------------------------------------------------
p_fig1 <- make_fig1_true_vs_observed(
  output_png = file.path(FIG_DIR, "Fig1_true_vs_observed.png"),
  output_rds = file.path(FIG_DIR, "Fig1_true_vs_observed.rds")
)
p_fig2 <- make_fig2_bayesian_estimation(
  output_png = file.path(FIG_DIR, "Fig2_bayesian_estimation.png"),
  output_rds = file.path(FIG_DIR, "Fig2_bayesian_estimation.rds")
)

#-------------------------------------------------------------------------------
# Figures 4 and 5: from the derived estimands
#-------------------------------------------------------------------------------
est_path = file.path(DERIVED_DIR, "estimands.rds")
if (!file.exists(est_path)) {
  stop("Missing ", est_path, ". Run scripts/03_derive_estimands.R first (Figures 1 and 2 were written).")
}
est <- readRDS(est_path)

p_trajectory <- compose_trajectory_figure(
  est$traj_data, est$obs_data,
  treat_results = est$treat_results,
  add_annotations = TRUE
)
ggplot2::ggsave(
  file.path(FIG_DIR, "Fig4_ogtt_trajectory_comparison.png"),
  plot = p_trajectory, width = 10, height = 5.5, dpi = 300
)
saveRDS(p_trajectory, file.path(FIG_DIR, "Fig4_ogtt_trajectory_comparison.rds"))

make_prior_sensitivity_forest(
  est$prior_sensitivity_df %>% dplyr::filter(prior_sd_scale != 999),
  output_path = file.path(FIG_DIR, "Fig5_prior_sensitivity_forest.png")
)

cat("Wrote Figures 1, 2, 4, 5 to", FIG_DIR, "\n")
