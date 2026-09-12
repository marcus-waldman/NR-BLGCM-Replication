#' Step 1: read the two OGTT visits, reshape to long bivariate format, add the
#' Legendre polynomial basis, join the visits, and restrict to participants
#' with a 36-week OGTT (N = 46). Writes derived/joined_01_long.rds, the
#' analysis data set used by every later step.
#'
#' Run from the repository root:  Rscript scripts/01_prepare_data.R
#' Requires the OGTT workbook (see data/raw/README.md); not distributed.

local({
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a) == 1) {
    setwd(normalizePath(file.path(dirname(sub("^--file=", "", a)), "..")))
  }
})
suppressMessages({
  library(tidyverse)
  library(orthopolynom)
})
source("R/_setup.R")

#-------------------------------------------------------------------------------
# Sheet 3 = baseline (31 weeks), suffix "0"; sheet 4 = post-intervention
# (36 weeks), suffix "1". Everyone with data on the sheet is read; the
# analytic sample is defined by the join below.
#-------------------------------------------------------------------------------
baseline_long <-
  read_ogtt_data(sheet = 3, path = OGTT_XLSX, require_completers = FALSE) %>%
  prepare_ogtt_long(suffix = "0") %>%
  add_legendre_basis()

post_long <-
  read_ogtt_data(sheet = 4, path = OGTT_XLSX, require_completers = FALSE) %>%
  prepare_ogtt_long(suffix = "1") %>%
  add_legendre_basis()

common_ptids <- intersect(baseline_long$PTID, post_long$PTID)
stopifnot(all(baseline_long$PTID %in% common_ptids))
stopifnot(all(post_long$PTID %in% common_ptids))

joined_01_long <- baseline_long |>
  dplyr::select(
    PTID, treat, minute_numeric, model_timepoint,
    t_std, dplyr::starts_with("leg_poly"),
    "logglucose0", "loginsulin0"
  ) %>%
  left_join(
    post_long %>% select(PTID, minute_numeric, "logglucose1", "loginsulin1"),
    by = c("PTID", "minute_numeric")
  )

# Keep participants with at least one 36-week glucose or insulin value.
PTIDs_keep <- joined_01_long %>%
  dplyr::group_by(PTID) %>%
  dplyr::summarise(all_miss = sum(is.na(logglucose1)) == 6 & sum(is.na(loginsulin1)) == 6) %>%
  dplyr::filter(!all_miss) %>%
  purrr::pluck("PTID")
if (length(PTIDs_keep) != N_EXPECTED) {
  stop("Unexpected sample size: ", length(PTIDs_keep), " (expected ", N_EXPECTED, ")")
}
joined_01_long <- joined_01_long %>% dplyr::filter(PTID %in% PTIDs_keep)

cat("Participants per arm:\n")
joined_01_long %>%
  dplyr::filter(minute_numeric == -15) %>%
  dplyr::count(treat) %>%
  print()

saveRDS(joined_01_long, file.path(DERIVED_DIR, "joined_01_long.rds"))
cat("Wrote", file.path(DERIVED_DIR, "joined_01_long.rds"),
    sprintf("(%d rows, %d columns)\n", nrow(joined_01_long), ncol(joined_01_long)))
