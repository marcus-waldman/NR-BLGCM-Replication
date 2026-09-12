#' Shared setup for every script in scripts/.
#'
#' Each script first sets the working directory to the repository root, then
#' sources this file. It defines the file locations, the literature constants
#' shared by the fitting and verification scripts, and sources every helper in
#' R/. Nothing here reads participant data.
#'
#' Environment overrides (so participant-level files never need to live inside
#' this repository):
#'   BLGCM_DATA_XLSX  path to the OGTT workbook   (default data/raw/OGTTdata.xlsx)
#'   BLGCM_FITS_DIR   directory holding the fits   (default fits/)

OGTT_XLSX = Sys.getenv("BLGCM_DATA_XLSX", unset = "data/raw/OGTTdata.xlsx")
FITS_DIR = Sys.getenv("BLGCM_FITS_DIR", unset = "fits")
DERIVED_DIR = "derived"
FIG_DIR = "figures"
dir.create(DERIVED_DIR, showWarnings = FALSE)
dir.create(FIG_DIR, showWarnings = FALSE)

# The five brms fits (primary + four sensitivity variants). Keys are used by
# scripts/02_fit_models.R, scripts/03_derive_estimands.R and
# scripts/verify_model_spec.R.
FIT_FILES <- c(
  primary = "fit_joined_01_v2.rds",
  tenthx = "fit_joined_tenthx_01.rds",
  diffuse_10x = "fit_diffuse_10x_01.rds",
  diffuse_100x = "fit_diffuse_100x_01.rds",
  indep_10x = "fit_indep_10x_01.rds"
)
fit_path <- function(key) {
  return(file.path(FITS_DIR, FIT_FILES[[key]]))
}

# Literature analytical CVs that center the residual-SD priors.
# Glucose: 0.58% (Aarsand et al., 2018, EuBIVAS).
# Insulin: 3.1%  (Carobene et al., 2022, EuBIVAS update).
GLUCOSE_CV_RESIDUAL_LIT = 0.0058
INSULIN_CV_RESIDUAL_LIT = 0.031

# Expected analytic sample: participants with both the 31-week and 36-week OGTT.
N_EXPECTED = 46

# Source every helper in R/ (this file excluded).
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  if (basename(f) != "_setup.R") {
    source(f)
  }
}
rm(f)
