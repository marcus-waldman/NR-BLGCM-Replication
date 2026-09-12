#' Run the whole pipeline in order. From the repository root:
#'   Rscript run_all.R
#'
#' Step 2 (model fitting) takes hours and is skipped for any fit already
#' present in FITS_DIR. Steps 3 and 5 take minutes.

local({
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a) == 1) {
    setwd(normalizePath(dirname(sub("^--file=", "", a))))
  }
})

steps <- c(
  "scripts/01_prepare_data.R",
  "scripts/02_fit_models.R",
  "scripts/03_derive_estimands.R",
  "scripts/04_make_figures.R",
  "scripts/05_make_figure6_priorsense.R"
)
rscript = file.path(R.home("bin"), "Rscript")
for (s in steps) {
  cat("\n################  ", s, "  ################\n")
  status = system2(rscript, shQuote(s))
  if (status != 0) {
    stop("Step failed: ", s, " (exit status ", status, ")")
  }
}
cat("\nPipeline complete. Figures are in figures/.\n")
