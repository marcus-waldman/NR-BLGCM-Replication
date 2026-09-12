#' Step 3: derive the quantities that Figures 4 and 5 plot from the fitted
#' models, and save them (fit-free) to derived/estimands.rds.
#'
#'   - AUC treatment effects (cross-adjusted ANCOVA, per posterior draw) for
#'     the primary fit in both modes (true-score, observed) and for the four
#'     sensitivity fits in true-score mode: six extractions.
#'   - prior_sensitivity_df   one row per analyte x prior-SD scale  (Figure 5)
#'   - traj_data / obs_data   36-week group trajectories               (Figure 4)
#'   - treat_results          AUC inset rows                            (Figure 4)
#'
#' Run from the repository root:  Rscript scripts/03_derive_estimands.R
#' Requires the five fits in FITS_DIR (step 2). The six extractions run on a
#' PSOCK cluster; each worker reads its fit from disk.

local({
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a) == 1) {
    setwd(normalizePath(file.path(dirname(sub("^--file=", "", a)), "..")))
  }
})
suppressMessages({
  library(tidyverse)
  library(brms)
  library(bayestestR)
  library(mice)
  library(pbapply)
})
source("R/_setup.R")

for (k in names(FIT_FILES)) {
  if (!file.exists(fit_path(k))) {
    stop("Missing fit '", k, "' at ", fit_path(k), ". Run scripts/02_fit_models.R first.")
  }
}

#-------------------------------------------------------------------------------
# 3.1  AUC treatment effects for all five fits
#-------------------------------------------------------------------------------
auc_te_jobs <- tibble::tribble(
  ~out_name, ~fit_key, ~mode,
  "auc_te_true", "primary", "true_score",
  "auc_te_obs", "primary", "observed",
  "auc_te_true_tenthx", "tenthx", "true_score",
  "auc_te_true_diffuse_10x", "diffuse_10x", "true_score",
  "auc_te_true_diffuse_100x", "diffuse_100x", "true_score",
  "auc_te_true_indep_10x", "indep_10x", "true_score"
)
auc_te_jobs$fit_path <- vapply(auc_te_jobs$fit_key, function(k) { return(normalizePath(fit_path(k))) }, character(1))

repo_root = getwd()
cl <- parallel::makeCluster(min(nrow(auc_te_jobs), max(1, parallel::detectCores() - 1)))
parallel::clusterExport(cl, "repo_root")
invisible(parallel::clusterEvalQ(cl, {
  setwd(repo_root)
  library(brms)
  library(dplyr)
  library(tidyr)
  source("R/legendre_polynomials.R")
  source("R/extract_clinical_indices.R")
}))
auc_te_results <- pbapply::pblapply(
  split(auc_te_jobs, seq_len(nrow(auc_te_jobs))),
  function(job) {
    fit <- readRDS(job$fit_path)
    return(extract_auc_treatment_effect(fit, mode = job$mode))
  },
  cl = cl
)
parallel::stopCluster(cl)
names(auc_te_results) <- auc_te_jobs$out_name

#-------------------------------------------------------------------------------
# 3.2  Prior-sensitivity table (Figure 5). prior_sd_scale 999 = independence
#      model; it is a structural alternative and is dropped from the figure.
#-------------------------------------------------------------------------------
auc_list2df <- function(in_list, prior_sd_scale) {
  out <- dplyr::bind_rows(
    with(in_list$AUC_glucose, data.frame(analyte = "Glucose", est = estimate, ci_lb = ci[1], ci_ub = ci[2])),
    with(in_list$AUC_insulin, data.frame(analyte = "Insulin", est = estimate, ci_lb = ci[1], ci_ub = ci[2]))
  )
  out$prior_sd_scale = prior_sd_scale
  return(out)
}
prior_sensitivity_df <- dplyr::bind_rows(
  auc_list2df(auc_te_results$auc_te_true_tenthx, 0.1),
  auc_list2df(auc_te_results$auc_te_true, 1),
  auc_list2df(auc_te_results$auc_te_true_diffuse_10x, 10),
  auc_list2df(auc_te_results$auc_te_true_diffuse_100x, 100),
  auc_list2df(auc_te_results$auc_te_true_indep_10x, 999)
) %>%
  dplyr::arrange(analyte, prior_sd_scale) %>%
  dplyr::mutate(ci_width = ci_ub - ci_lb) %>%
  dplyr::relocate(analyte, prior_sd_scale)

#-------------------------------------------------------------------------------
# 3.3  36-week trajectories + AUC inset rows (Figure 4), primary fit
#-------------------------------------------------------------------------------
fit_joined_01 <- readRDS(fit_path("primary"))
traj_data <- compute_true_score_trajectories(fit_joined_01, suffixes = "1")
obs_data <- compute_observed_group_means(fit_joined_01, suffixes = "1")

# The "True-score" label is required by plot_treatment_inset(), which greps
# on it to split true-score from observed rows.
make_treat_row <- function(index, method, est_ci) {
  return(data.frame(
    index = index,
    method = method,
    estimate = est_ci$estimate,
    conf.low = est_ci$ci[1],
    conf.high = est_ci$ci[2],
    stringsAsFactors = FALSE
  ))
}
treat_results <- dplyr::bind_rows(
  make_treat_row("Glucose AUC", "Observed\n(Rubin's rules)", auc_te_results$auc_te_obs$AUC_glucose),
  make_treat_row("Glucose AUC", "True-score\n(ME-corrected)", auc_te_results$auc_te_true$AUC_glucose),
  make_treat_row("Insulin AUC", "Observed\n(Rubin's rules)", auc_te_results$auc_te_obs$AUC_insulin),
  make_treat_row("Insulin AUC", "True-score\n(ME-corrected)", auc_te_results$auc_te_true$AUC_insulin)
)

#-------------------------------------------------------------------------------
# 3.4  Save + console report
#-------------------------------------------------------------------------------
estimands <- list(
  auc_te = auc_te_results,
  prior_sensitivity_df = prior_sensitivity_df,
  traj_data = traj_data,
  obs_data = obs_data,
  treat_results = treat_results
)
saveRDS(estimands, file.path(DERIVED_DIR, "estimands.rds"))

cat("\n== AUC treatment effects (CHOICE - LC/Conventional), primary fit ==\n")
print(treat_results, row.names = FALSE)
cat("\n== Prior sensitivity (true-score mode) ==\n")
print(prior_sensitivity_df, row.names = FALSE)
cat("\nWrote", file.path(DERIVED_DIR, "estimands.rds"), "\n")
