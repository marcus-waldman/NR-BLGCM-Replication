#' Step 2: fit the primary B-LGCM and the four sensitivity variants with brms
#' (cmdstanr backend). Each fit is skipped when its .rds already exists in
#' FITS_DIR, so a reader who has the fits can go straight to step 3.
#'
#' Run from the repository root:  Rscript scripts/02_fit_models.R
#' Runtime: hours per fit (4 chains x 3000 iterations, adapt_delta 0.99,
#' max_treedepth 25, within-chain threading). Fits are ~135 MB each and hold
#' participant-level rows in fit$data, so keep them out of version control.
#'
#' Fits (file names in R/_setup.R):
#'   primary       informative priors, prior_sd_scale = 1
#'   tenthx        prior SDs scaled 0.1x
#'   diffuse_10x   prior SDs scaled 10x
#'   diffuse_100x  prior SDs scaled 100x
#'   indep_10x     independence model (no cross-response random-effect
#'                 correlations), prior SDs scaled 10x

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
dir.create(FITS_DIR, showWarnings = FALSE, recursive = TRUE)

joined_01_long <- readRDS(file.path(DERIVED_DIR, "joined_01_long.rds"))

# Sampler settings shared by all five fits.
SEED = 42
CHAINS = 4
WARMUP = 1000
ITER = 3000
CONTROL <- list(adapt_delta = 0.99, max_treedepth = 25)

#-------------------------------------------------------------------------------
# Model specification: formula + priors (see R/likelihood_specification.R and
# R/prior_specification.R). The Legendre-coefficient priors are derived at
# runtime from the literature CSVs in priors/.
#-------------------------------------------------------------------------------
mvform_joined_01 <- create_joined_growth_curve_formula_n(suffixes = c("0", "1"), rescor = FALSE)
print(mvform_joined_01)

make_priors <- function(prior_sd_scale = 1, independence = FALSE) {
  return(create_joined_growth_curve_priors_n(
    suffixes = c("0", "1"),
    rescor = FALSE,
    independence = independence,
    data = joined_01_long,
    prior_sd_scale = prior_sd_scale,
    glucose_cv_residual = GLUCOSE_CV_RESIDUAL_LIT,
    insulin_cv_residual = INSULIN_CV_RESIDUAL_LIT
  ))
}
dpriors_joined_01 <- make_priors(prior_sd_scale = 1)
print(dpriors_joined_01)

#-------------------------------------------------------------------------------
# 2.1  Primary fit
#-------------------------------------------------------------------------------
if (!file.exists(fit_path("primary"))) {
  fit_joined_01 <- brms::brm(
    formula = mvform_joined_01,
    data = joined_01_long,
    prior = dpriors_joined_01,
    backend = "cmdstanr",
    control = CONTROL,
    seed = SEED,
    chains = CHAINS,
    cores = CHAINS,
    threads = threading(4),
    warmup = WARMUP,
    iter = ITER,
    init = 1,
    file = sub("\\.rds$", "", fit_path("primary")),
    file_compress = TRUE,
    file_refit = "always"
  )
} else {
  cat("Primary fit exists, loading:", fit_path("primary"), "\n")
  fit_joined_01 <- readRDS(fit_path("primary"))
}

#-------------------------------------------------------------------------------
# 2.2  Per-chain inits for the prior-scaled refits, drawn from the primary
#      posterior (the parameter block is identical across these refits).
#-------------------------------------------------------------------------------
sens_inits_joined <- make_inits_from_brms(fit_joined_01, n_chains = CHAINS)

fit_sensitivity <- function(key, prior_sd_scale, threads_per_chain) {
  path = fit_path(key)
  if (file.exists(path)) {
    cat("Fit exists, skipping:", path, "\n")
    return(invisible(NULL))
  }
  fit <- brms::brm(
    formula = mvform_joined_01,
    data = joined_01_long,
    prior = make_priors(prior_sd_scale = prior_sd_scale),
    backend = "cmdstanr",
    control = CONTROL,
    seed = SEED,
    chains = CHAINS,
    cores = CHAINS,
    threads = threading(threads_per_chain),
    warmup = WARMUP,
    iter = ITER,
    init = sens_inits_joined,
    file = sub("\\.rds$", "", path),
    file_compress = TRUE,
    file_refit = "always"
  )
  return(invisible(fit))
}

#-------------------------------------------------------------------------------
# 2.3  Prior-sensitivity refits: 0.1x, 10x, 100x prior SDs
#-------------------------------------------------------------------------------
fit_sensitivity("tenthx", prior_sd_scale = 0.1, threads_per_chain = 2)
fit_sensitivity("diffuse_10x", prior_sd_scale = 10, threads_per_chain = 4)
fit_sensitivity("diffuse_100x", prior_sd_scale = 100, threads_per_chain = 4)

#-------------------------------------------------------------------------------
# 2.4  Structural alternative: independence model, 10x prior SDs
#-------------------------------------------------------------------------------
if (!file.exists(fit_path("indep_10x"))) {
  mvform_indep_01 <- create_joined_growth_curve_formula_n(
    suffixes = c("0", "1"), rescor = FALSE, independence = TRUE
  )
  fit_indep_10x_01 <- brms::brm(
    formula = mvform_indep_01,
    data = joined_01_long,
    prior = make_priors(prior_sd_scale = 10, independence = TRUE),
    backend = "cmdstanr",
    control = CONTROL,
    seed = SEED,
    chains = CHAINS,
    cores = CHAINS,
    threads = threading(4),
    warmup = WARMUP,
    iter = ITER,
    init = 1 / 10,
    file = sub("\\.rds$", "", fit_path("indep_10x")),
    file_compress = TRUE,
    file_refit = "always"
  )
} else {
  cat("Fit exists, skipping:", fit_path("indep_10x"), "\n")
}

cat("All five fits present in", FITS_DIR, "\n")
