#' Transform Cleaned Data to Long Format with Bivariate Structure
#'
#' Transforms raw OGTT data to long format suitable for growth curve modeling:
#'   1. Selects and renames insulin and glucose columns
#'   2. Pivots to long format
#'   3. Maps minute measurements to model timepoints (-15 and 0 → 0)
#'   4. Joins insulin and glucose into bivariate structure
#'   5. Log transforms both analytes
#'   6. Returns log-transformed bivariate data
#'
#' @param dat Cleaned data from read_ogtt_data()
#' @param suffix Character suffix for variable names ("0" = baseline, "1" = post, "2" = postpartum)
#'
#' @return Long format data with columns:
#'   - PTID, treat, minute_numeric, model_timepoint
#'   - glucoseX, insulinX (raw values)
#'   - logglucoseX, loginsulinX (log-transformed)
#'
#' @examples
#' baseline_long <- read_ogtt_data(3, proj_wd) |> prepare_ogtt_long(suffix = "0")
#'
#' @export
prepare_ogtt_long <- function(dat, suffix = "0") {
  # Dynamic column names
  glu_col <- paste0("glucose", suffix)
  ins_col <- paste0("insulin", suffix)
  logglu_col <- paste0("logglucose", suffix)
  logins_col <- paste0("loginsulin", suffix)
  
  # Select and rename columns
  dat_selected <- dat |>
    select(
      PTID, treat,
      INS00 = `Insulin -15`, INS0 = `Fasting Insulin`,
      INS30 = `Insulin30`, INS60 = `Insulin60`, INS90 = `Insulin90`, INS120 = `Insulin120`,
      GLU00 = `Glucose-15`, GLU0 = `Fasting Glucose`,
      GLU30 = `Glucose30`, GLU60 = `Glucose60`, GLU90 = `Glucose90`, GLU120 = `Glucose120`
    )
  
  print("dat_selected:")
  head(dat_selected)
  # Pivot to long format
  dat_long <- dat_selected |>
    pivot_longer(INS00:GLU120, names_to = "name", values_to = "y") |>
    mutate(
      outcome = ifelse(startsWith(name, "INS"), "insulin", "glucose"),
      minute_str = gsub("INS|GLU", "", name),
      minute_numeric = case_when(
        minute_str == "00" ~ -15,
        minute_str == "0" ~ 0,
        minute_str == "30" ~ 30,
        minute_str == "60" ~ 60,
        minute_str == "90" ~ 90,
        minute_str == "120" ~ 120,
        TRUE ~ NA_real_
      ),
    # Map strings to numeric with -15 timepoint -> t = 0
      model_timepoint = case_when(
        minute_numeric == -15 ~ 0,
        minute_numeric == 0 ~ 0,
        minute_numeric == 30 ~ 30,
        minute_numeric == 60 ~ 60,
        minute_numeric == 90 ~ 90,
        minute_numeric == 120 ~ 120,
        TRUE ~ NA_real_
      )
    )
  print("dat_long")
  head(dat_long)
  
  
  # Join insulin and glucose, rename with suffix
  dat_joined <- dat_long |>
    filter(outcome == "insulin") |>
    select(PTID, treat, !!ins_col := y, minute_numeric, model_timepoint) |>
    left_join(
      dat_long |> filter(outcome == "glucose") |> select(PTID, minute_numeric, !!glu_col := y),
      by = c("PTID", "minute_numeric")
    ) |>
    arrange(PTID, minute_numeric)
  
  # Log transform
  dat_joined <- dat_joined |>
    mutate(
      !!logglu_col := log(.data[[glu_col]]),
      !!logins_col := log(.data[[ins_col]])
    )
  print("dat_joined")
  print(dat_joined)
  
  return(dat_joined)
}
