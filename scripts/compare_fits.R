#' Refit check: does a fresh run of scripts/02_fit_models.R reproduce the
#' reference fits?
#'
#' The fits use within-chain threading with dynamic scheduling, so posterior
#' draws are not bitwise reproducible even with the same seed. The check is
#' therefore statistical: for every fit present in both directories,
#'   1. per-parameter z = (mean_new - mean_ref) / sqrt(mcse_ref^2 + mcse_new^2)
#'      across all monitored parameters; report max |z| and the share above 3.
#'      Under exact reproduction up to Monte Carlo error, |z| behaves like a
#'      standard normal, so with ~2,000 parameters a max around 3.5-4 is
#'      expected and anything above Z_FAIL flags a real difference.
#'   2. posterior-SD ratio new/ref, min and max.
#'   3. for the primary fit, the AUC treatment effects (true-score and
#'      observed modes) side by side, with the difference in estimates as a
#'      fraction of the reference interval width.
#'
#' Run from the repository root:
#'   Rscript scripts/compare_fits.R --ref <reference fits dir> --new <fresh fits dir> [--keys primary,tenthx,...]
#' Exits with status 1 if any fit fails the z or SD-ratio thresholds.

local({
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a) == 1) {
    setwd(normalizePath(file.path(dirname(sub("^--file=", "", a)), "..")))
  }
})
suppressMessages({
  library(tidyverse)
  library(brms)
  library(posterior)
  library(bayestestR)
  library(mice)
})
source("R/_setup.R")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i = match(flag, args)
  if (is.na(i) || i == length(args)) {
    return(default)
  }
  return(args[i + 1])
}
REF_DIR = get_arg("--ref")
NEW_DIR = get_arg("--new")
KEYS <- strsplit(get_arg("--keys", paste(names(FIT_FILES), collapse = ",")), ",")[[1]]
if (is.null(REF_DIR) || is.null(NEW_DIR)) {
  stop("Usage: Rscript scripts/compare_fits.R --ref <dir> --new <dir> [--keys k1,k2]")
}
Z_FAIL = 5
SD_RATIO_FAIL <- c(0.8, 1.25)
AUC_FRAC_FAIL = 0.10

summarise_fit <- function(fit) {
  d <- posterior::as_draws_array(fit)
  s <- posterior::summarise_draws(d, "mean", "sd", posterior::default_mcse_measures())
  return(s[, c("variable", "mean", "sd", "mcse_mean")])
}

n_fail = 0
for (key in KEYS) {
  p_ref = file.path(REF_DIR, FIT_FILES[[key]])
  p_new = file.path(NEW_DIR, FIT_FILES[[key]])
  cat(sprintf("\n==== %s ====\n  ref: %s\n  new: %s\n", key, p_ref, p_new))
  if (!file.exists(p_ref) || !file.exists(p_new)) {
    cat("  SKIP: missing on one side\n")
    next
  }
  fit_ref <- readRDS(p_ref)
  fit_new <- readRDS(p_new)

  s_ref <- summarise_fit(fit_ref)
  s_new <- summarise_fit(fit_new)
  both <- dplyr::inner_join(s_ref, s_new, by = "variable", suffix = c("_ref", "_new")) %>%
    dplyr::filter(is.finite(sd_ref), is.finite(sd_new), sd_ref > 0, sd_new > 0,
                  !grepl("^lp__$|^lprior$", variable)) %>%
    dplyr::mutate(
      z = (mean_new - mean_ref) / sqrt(mcse_mean_ref^2 + mcse_mean_new^2),
      sd_ratio = sd_new / sd_ref
    )
  n_par = nrow(both)
  max_z = max(abs(both$z))
  share_z3 = mean(abs(both$z) > 3)
  sd_rng <- range(both$sd_ratio)
  cat(sprintf("  parameters compared     : %d (of %d ref / %d new)\n", n_par, nrow(s_ref), nrow(s_new)))
  cat(sprintf("  max |z| on posterior mean: %.2f   (share |z|>3: %.4f)\n", max_z, share_z3))
  cat(sprintf("  posterior-SD ratio new/ref: [%.3f, %.3f]\n", sd_rng[1], sd_rng[2]))
  worst <- both %>% dplyr::arrange(dplyr::desc(abs(z))) %>% dplyr::slice_head(n = 5)
  cat("  largest |z|:\n")
  print(as.data.frame(worst[, c("variable", "mean_ref", "mean_new", "z", "sd_ratio")]), row.names = FALSE, digits = 4)

  ok_z = max_z < Z_FAIL
  ok_sd = sd_rng[1] > SD_RATIO_FAIL[1] && sd_rng[2] < SD_RATIO_FAIL[2]
  cat(sprintf("  z check   : %s (threshold %g)\n", if (ok_z) { "PASS" } else { "FAIL" }, Z_FAIL))
  cat(sprintf("  SD check  : %s (window %g-%g)\n", if (ok_sd) { "PASS" } else { "FAIL" }, SD_RATIO_FAIL[1], SD_RATIO_FAIL[2]))

  ok_auc = TRUE
  if (key == "primary") {
    cat("\n  AUC treatment effects (CHOICE - LC/Conventional), primary fit:\n")
    rows <- list()
    for (mode in c("true_score", "observed")) {
      r_ref <- extract_auc_treatment_effect(fit_ref, mode = mode)
      r_new <- extract_auc_treatment_effect(fit_new, mode = mode)
      for (an in c("AUC_glucose", "AUC_insulin")) {
        w_ref = diff(r_ref[[an]]$ci)
        rows[[length(rows) + 1]] <- data.frame(
          mode = mode, analyte = an,
          est_ref = r_ref[[an]]$estimate, est_new = r_new[[an]]$estimate,
          lo_ref = r_ref[[an]]$ci[1], lo_new = r_new[[an]]$ci[1],
          hi_ref = r_ref[[an]]$ci[2], hi_new = r_new[[an]]$ci[2],
          d_est_frac = abs(r_new[[an]]$estimate - r_ref[[an]]$estimate) / w_ref,
          d_width_frac = abs(diff(r_new[[an]]$ci) - w_ref) / w_ref
        )
      }
    }
    auc_tab <- dplyr::bind_rows(rows)
    print(auc_tab, row.names = FALSE, digits = 4)
    ok_auc = all(auc_tab$d_est_frac < AUC_FRAC_FAIL) && all(auc_tab$d_width_frac < AUC_FRAC_FAIL)
    cat(sprintf("  AUC check : %s (estimate and width shifts < %g of reference CrI width)\n",
                if (ok_auc) { "PASS" } else { "FAIL" }, AUC_FRAC_FAIL))
  }

  if (!(ok_z && ok_sd && ok_auc)) {
    n_fail = n_fail + 1
  }
  rm(fit_ref, fit_new)
  invisible(gc())
}

cat(sprintf("\n%s\n", if (n_fail == 0) { "ALL COMPARED FITS REPRODUCE WITHIN MONTE CARLO ERROR" } else { paste(n_fail, "FIT(S) FAILED") }))
if (n_fail > 0) {
  quit(status = 1)
}
