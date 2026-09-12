#' Read and Clean Raw OGTT Data from Excel
#'
#' Reads OGTT data from specified Excel sheet and performs initial cleaning:
#'   - Optionally filters to completers only
#'   - Converts PTID to factor
#'   - Converts Group to factor with labels
#'   - Creates binary treatment indicator
#'
#' @param sheet Integer sheet number (2 = postpartum, 3 = baseline, 4 = post-intervention)
#' @param path Character path to the OGTT Excel workbook (see data/raw/README.md
#'   for the expected layout). Defaults to `OGTT_XLSX`, set in R/_setup.R.
#' @param require_completers Logical, filter to study completers (default: TRUE).
#'   Set to FALSE to include all participants with data on this sheet,
#'   regardless of whether they completed all study phases.
#'
#' @return Data frame with cleaned OGTT data
#'
#' @examples
#' # Default: only 3-phase completers (n=44)
#' baseline_raw <- read_ogtt_data(sheet = 3, path = "data/raw/OGTTdata.xlsx")
#'
#' # All participants with baseline data
#' baseline_all <- read_ogtt_data(sheet = 3, path = "data/raw/OGTTdata.xlsx",
#'                                require_completers = FALSE)
#'
#' @export
read_ogtt_data <- function(sheet, path = OGTT_XLSX, require_completers = TRUE) {
  if (!file.exists(path)) {
    stop("OGTT workbook not found at '", path, "'. The participant-level data are ",
         "not distributed with this repository; see data/raw/README.md.")
  }
  dat <- readxl::read_excel(path, sheet = sheet)

  if (require_completers) {
    dat <- dat |> filter(!is.na(Completers))
  }

  dat = dat %>%
    mutate(
      PTID = as.factor(PTID),
      Group = factor(Group, levels = c(0, 1), labels = c("LC/Conventional", "CHOICE")),
      treat = as.integer(Group == "CHOICE")
    )

  return(dat)
}
