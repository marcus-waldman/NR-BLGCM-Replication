#' Verification: the model specification this repository rebuilds matches the
#' saved fits exactly.
#'
#' For each of the five fits in FITS_DIR, rebuild the brms formula, priors and
#' data from scripts 01/02 and compare against what the fit object stores:
#'   1. Stan program text   brms::make_stancode(...)  vs  brms::stancode(fit)
#'      (prior values are embedded in the program, so this ties the priors)
#'   2. Stan data list      brms::make_standata(...)  vs  brms::standata(fit)
#'      (ties the response vectors, design matrices, missingness indices)
#'   3. Model frame         columns of fit$data       vs  derived/joined_01_long.rds
#' `grainsize` is reported separately: it depends only on the threading
#' setting and does not affect the posterior.
#'
#' Run from the repository root:  Rscript scripts/verify_model_spec.R
#' Exits with status 1 on any mismatch.

local({
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a) == 1) {
    setwd(normalizePath(file.path(dirname(sub("^--file=", "", a)), "..")))
  }
})
suppressMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
})
source("R/_setup.R")

joined_01_long <- readRDS(file.path(DERIVED_DIR, "joined_01_long.rds"))

specs <- list(
  primary = list(prior_sd_scale = 1, independence = FALSE, threads = 4),
  tenthx = list(prior_sd_scale = 0.1, independence = FALSE, threads = 2),
  diffuse_10x = list(prior_sd_scale = 10, independence = FALSE, threads = 4),
  diffuse_100x = list(prior_sd_scale = 100, independence = FALSE, threads = 4),
  indep_10x = list(prior_sd_scale = 10, independence = TRUE, threads = 4)
)

compare_standata <- function(new, old) {
  keys <- union(names(new), names(old))
  diffs <- character(0)
  for (k in keys) {
    if (!(k %in% names(new)) || !(k %in% names(old))) {
      diffs <- c(diffs, paste0(k, " (missing on one side)"))
      next
    }
    same = isTRUE(all.equal(unname(new[[k]]), unname(old[[k]]), check.attributes = FALSE))
    if (!same) {
      diffs <- c(diffs, k)
    }
  }
  return(diffs)
}

n_fail = 0
for (key in names(specs)) {
  sp <- specs[[key]]
  cat(sprintf("\n==== %s (%s) ====\n", key, fit_path(key)))
  if (!file.exists(fit_path(key))) {
    cat("  SKIP: fit not found\n")
    next
  }
  fit <- readRDS(fit_path(key))

  form <- create_joined_growth_curve_formula_n(
    suffixes = c("0", "1"), rescor = FALSE, independence = sp$independence
  )
  pri <- create_joined_growth_curve_priors_n(
    suffixes = c("0", "1"), rescor = FALSE, independence = sp$independence,
    data = joined_01_long, prior_sd_scale = sp$prior_sd_scale,
    glucose_cv_residual = GLUCOSE_CV_RESIDUAL_LIT,
    insulin_cv_residual = INSULIN_CV_RESIDUAL_LIT
  )

  # 1. Stan program
  code_new = brms::make_stancode(form, data = joined_01_long, prior = pri,
                                 backend = "cmdstanr", threads = threading(sp$threads))
  code_old = brms::stancode(fit)
  ok_code = identical(as.character(code_new), as.character(code_old))
  cat(sprintf("  Stan program identical : %s\n", ok_code))
  if (!ok_code) {
    ln_new <- strsplit(as.character(code_new), "\n")[[1]]
    ln_old <- strsplit(as.character(code_old), "\n")[[1]]
    first = which(ln_new != ln_old)[1]
    cat("    first differing line", first, "\n    new:", ln_new[first], "\n    old:", ln_old[first], "\n")
  }

  # 2. Stan data
  sd_new <- brms::make_standata(form, data = joined_01_long, prior = pri,
                                backend = "cmdstanr", threads = threading(sp$threads))
  sd_old <- brms::standata(fit)
  diffs <- compare_standata(sd_new, sd_old)
  grain_diff = "grainsize" %in% diffs
  diffs <- setdiff(diffs, "grainsize")
  ok_data = length(diffs) == 0
  cat(sprintf("  Stan data identical    : %s%s\n", ok_data,
              if (grain_diff) { "  (grainsize differs: threading-only, posterior-irrelevant)" } else { "" }))
  if (!ok_data) {
    cat("    differing elements:", paste(diffs, collapse = ", "), "\n")
  }

  # 3. Model frame: every column brms kept must match the rebuilt data
  shared <- intersect(names(fit$data), names(joined_01_long))
  ok_frame = isTRUE(all.equal(
    as.data.frame(fit$data)[, shared],
    as.data.frame(joined_01_long)[, shared],
    check.attributes = FALSE
  )) && length(shared) == ncol(fit$data)
  cat(sprintf("  Model frame identical  : %s  (%d shared columns, %d rows)\n",
              ok_frame, length(shared), nrow(fit$data)))

  if (!(ok_code && ok_data && ok_frame)) {
    n_fail = n_fail + 1
  }
  rm(fit)
  invisible(gc())
}

cat(sprintf("\n%s\n", if (n_fail == 0) { "ALL FITS MATCH THE REBUILT SPECIFICATION" } else { paste(n_fail, "FIT(S) DO NOT MATCH") }))
if (n_fail > 0) {
  quit(status = 1)
}
