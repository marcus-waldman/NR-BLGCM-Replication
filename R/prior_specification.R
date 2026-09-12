#' Default Paths to External OGTT Reference CSVs (v2)
#'
#' The v2 CSVs preserve the two independent visual-rater readings of the
#' published figures, allowing rater noise to enter the prior via the lm
#' standard error. Pregnancy is the only biological state used in this
#' repository (baseline 31wk + post 36wk); postpartum priors have been removed.
#'
#' @keywords internal
external_ogtt_csv_paths <- function() {
  list(
    glucose = file.path(
      "priors", "carlson-log-glucose",
      "external_ogtt_glucose_reference v2.csv"
    ),
    insulin = file.path(
      "priors", "mittendorfer-log-insulin",
      "mittendorfer_ogtt_insulin_reference v2.csv"
    )
  )
}

#' Load External OGTT Reference (v2) and Reshape to Long Per-Rater Format
#'
#' Handles both v2 schemas:
#'   - Carlson (glucose): wide format with columns
#'     `gdm_my_estimate`, `gdm_other_estimate`, `no_gdm_my_estimate`,
#'     `no_gdm_other_estimate` plus a precomputed `*_avg` (ignored here).
#'   - Mittendorfer (insulin): long format with columns `time_min`, `rater`,
#'     `no_gdm_35wk`, `gdm_35wk`.
#'
#' @param analyte "glucose" or "insulin".
#' @param csv_path Optional override for the CSV path. When NULL (default), uses
#'   `external_ogtt_csv_paths()` resolved against the current working directory.
#'
#' @return Data frame with columns: `time_min`, `rater`, `group` ("gdm" or
#'   "no_gdm"), `value` (raw analyte concentration). One row per
#'   timepoint × rater × group.
#'
#' @keywords internal
load_external_ogtt_long <- function(analyte, csv_path = NULL) {
  if (is.null(csv_path)) {
    csv_path <- external_ogtt_csv_paths()[[analyte]]
  }
  if (!file.exists(csv_path)) {
    stop("External OGTT CSV not found: ", csv_path)
  }
  raw <- utils::read.csv(csv_path, check.names = TRUE, stringsAsFactors = FALSE)

  if (analyte == "glucose") {
    # Carlson v2: wide rater columns. Rename to a uniform per-rater long table.
    needed <- c("time_min",
                "gdm_my_estimate", "gdm_other_estimate",
                "no_gdm_my_estimate", "no_gdm_other_estimate")
    missing <- setdiff(needed, names(raw))
    if (length(missing) > 0) {
      stop("Glucose v2 CSV is missing required columns: ",
           paste(missing, collapse = ", "))
    }
    long <- rbind(
      data.frame(time_min = raw$time_min, rater = "rater_1",
                 group = "gdm",    value = raw$gdm_my_estimate),
      data.frame(time_min = raw$time_min, rater = "rater_2",
                 group = "gdm",    value = raw$gdm_other_estimate),
      data.frame(time_min = raw$time_min, rater = "rater_1",
                 group = "no_gdm", value = raw$no_gdm_my_estimate),
      data.frame(time_min = raw$time_min, rater = "rater_2",
                 group = "no_gdm", value = raw$no_gdm_other_estimate)
    )
  } else if (analyte == "insulin") {
    # Mittendorfer v2: already long over rater; pivot the two group columns.
    needed <- c("time_min", "rater", "no_gdm_35wk", "gdm_35wk")
    missing <- setdiff(needed, names(raw))
    if (length(missing) > 0) {
      stop("Insulin v2 CSV is missing required columns: ",
           paste(missing, collapse = ", "))
    }
    long <- rbind(
      data.frame(time_min = raw$time_min, rater = raw$rater,
                 group = "gdm",    value = raw$gdm_35wk),
      data.frame(time_min = raw$time_min, rater = raw$rater,
                 group = "no_gdm", value = raw$no_gdm_35wk)
    )
  } else {
    stop("analyte must be 'glucose' or 'insulin'")
  }

  # Restrict to OGTT window [0, 120] min that maps to Legendre domain [-1, 1].
  long <- long[long$time_min >= 0 & long$time_min <= 120, , drop = FALSE]
  long$t_scaled <- long$time_min / 60 - 1
  long$log_value <- log(long$value)
  long
}

#' Standard Legendre Polynomial Basis P0..P4
#' @keywords internal
legendre_basis_p0_p4 <- function(x) {
  data.frame(
    P0 = rep(1, length(x)),
    P1 = x,
    P2 = (3 * x^2 - 1) / 2,
    P3 = (5 * x^3 - 3 * x) / 2,
    P4 = (35 * x^4 - 30 * x^2 + 3) / 8
  )
}

#' Derive Empirical Legendre Priors from a v2 External-OGTT CSV
#'
#' For the requested `analyte`, this:
#'   1. Loads the v2 CSV with both raters preserved.
#'   2. For each external group (GDM, No-GDM), fits a Legendre P0..P4 lm to
#'      log(value) using stacked rater × timepoint rows. The lm standard error
#'      on each coefficient therefore reflects digitization / rater noise.
#'   3. Centers the prior on the GDM coefficient and constructs the prior SD by
#'      combining rater noise (SE_gdm) and between-group spread
#'      (|coef_gdm - coef_no_gdm| / 1.5) in quadrature.
#'
#' @param analyte "glucose" or "insulin".
#' @param csv_path Optional override for the v2 CSV path.
#'
#' @return Data frame with one row per term (P0..P4), columns:
#'   `term`, `coef_gdm`, `se_gdm`, `coef_no_gdm`, `se_no_gdm`,
#'   `group_spread`, `prior_mean`, `prior_sd`.
#'
#' @keywords internal
derive_legendre_priors_from_csv <- function(analyte, csv_path = NULL) {
  long <- load_external_ogtt_long(analyte, csv_path = csv_path)
  basis <- legendre_basis_p0_p4(long$t_scaled)
  dat <- cbind(long, basis)

  fit_group <- function(g) {
    sub <- dat[dat$group == g, , drop = FALSE]
    fit <- stats::lm(log_value ~ P0 + P1 + P2 + P3 + P4 - 1, data = sub)
    smry <- summary(fit)$coefficients
    data.frame(
      term = rownames(smry),
      coef = smry[, "Estimate"],
      se   = smry[, "Std. Error"],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }

  gdm    <- fit_group("gdm")
  no_gdm <- fit_group("no_gdm")
  names(gdm)    <- c("term", "coef_gdm",    "se_gdm")
  names(no_gdm) <- c("term", "coef_no_gdm", "se_no_gdm")

  out <- merge(gdm, no_gdm, by = "term", sort = FALSE)
  out$group_spread <- abs(out$coef_gdm - out$coef_no_gdm) / 1.5
  out$prior_mean   <- out$coef_gdm
  out$prior_sd     <- sqrt(out$se_gdm^2 + out$group_spread^2)
  out
}

#' Create Empirical Priors for Legendre Polynomial Coefficients
#'
#' Generates priors for Legendre polynomial coefficients (P1-P4) by fitting
#' Legendre polynomials at runtime to the v2 external OGTT reference CSVs.
#' Pregnancy is the only biological state supported in this repository.
#'
#' @param suffix Character suffix for timepoint ("0" = baseline 31wk,
#'   "1" = post-intervention 36wk). Used only to name the brms response
#'   variable; the external prior values themselves do not vary by suffix.
#' @param analyte Character, "glucose" or "insulin".
#' @param scale Positive scalar multiplier applied to every prior SD (default
#'   1.0). Means are unchanged; only the SD of each Normal prior is multiplied
#'   by `scale`. Use values > 1 to widen the empirical priors (e.g.,
#'   `scale = 10` for diffuse sensitivity refits).
#' @param csv_path Optional override for the v2 CSV path; defaults to the
#'   appropriate file under `priors/` (see
#'   `external_ogtt_csv_paths()`).
#'
#' @return List of brms priors for P1-P4 coefficients (control and treatment
#'   groups).
#'
#' @details
#' Empirical prior sources (pregnancy only):
#'   - Glucose: Carlson 2025 (`external_ogtt_glucose_reference v2.csv`,
#'     two raters preserved as separate columns; GDM-centered).
#'   - Insulin: Mittendorfer 2023
#'     (`mittendorfer_ogtt_insulin_reference v2.csv`, two raters in long
#'     format, 35-week groups; GDM-centered).
#'
#' Prior construction per coefficient P1..P4:
#'   - For each group (GDM, No-GDM), fit `lm(log(value) ~ P0+P1+P2+P3+P4 - 1)`
#'     with stacked rater × timepoint rows. The lm standard error therefore
#'     reflects digitization / rater noise.
#'   - `prior_mean = coef_gdm`
#'   - `prior_sd = sqrt(se_gdm^2 + ((coef_gdm - coef_no_gdm)/1.5)^2)`
#'     (rater noise + between-group spread combined in quadrature).
#'
#' P0 (intercept) receives a data-scaled prior elsewhere in
#' `create_joined_growth_curve_priors_n()` and is not handled here. Both
#' control and treatment groups receive identical empirical priors since
#' external OGTT data does not distinguish treatment arms.
#'
#' @examples
#' # Pregnancy priors for glucose at baseline (suffix "0")
#' priors_glucose_pre <- create_legendre_coefficient_priors("0", "glucose")
#'
#' @export
create_legendre_coefficient_priors <- function(suffix, analyte,
                                               scale = 1.0,
                                               csv_path = NULL) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' required but not installed.")
  }
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("scale must be a positive finite scalar.")
  }
  if (!(analyte %in% c("glucose", "insulin"))) {
    stop("analyte must be 'glucose' or 'insulin'")
  }

  outcome_var <- paste0("log", analyte, suffix)
  prior_tbl <- derive_legendre_priors_from_csv(analyte, csv_path = csv_path)

  # P0 is data-scaled in create_joined_growth_curve_priors_n(); skip here.
  terms_to_emit <- c("P1", "P2", "P3", "P4")
  fmt <- function(x) signif(x, 4)

  priors <- list()
  for (tm in terms_to_emit) {
    row <- prior_tbl[prior_tbl$term == tm, , drop = FALSE]
    if (nrow(row) != 1L) {
      stop("Could not locate Legendre term '", tm, "' for analyte '", analyte,
           "' in derived prior table.")
    }
    mu <- fmt(row$prior_mean)
    sd <- fmt(row$prior_sd * scale)
    coef_idx <- sub("^P", "", tm)  # "1", "2", "3", "4"
    for (grp in c("ctrl", "treat")) {
      priors <- c(priors, list(brms::set_prior(
        sprintf("normal(%s, %s)", mu, sd),
        class = "b",
        coef  = sprintf("leg_poly_%s_%s", coef_idx, grp),
        resp  = outcome_var
      )))
    }
  }
  do.call(c, priors)
}

#' Create Priors for N-Timepoint Joined Growth Curve Model
#'
#' Generalized version of create_joined_growth_curve_priors() that handles
#' arbitrary number of timepoints. Generates priors for 2N outcomes
#' (N timepoints × 2 analytes).
#'
#' @param suffixes Character vector of timepoint suffixes (e.g., c("0", "1", "2"))
#' @param what Character vector specifying which prior parameter classes to include.
#'   Options: "beta" (fixed effects), "sigma" (residual SD intercepts), "wiggle"
#'   (residual SD slopes), "rescor" (residual correlations), "sd" (random effect SDs),
#'   "etacor" (random effect correlations). Default: all classes.
#' @param rescor Logical, whether the formula uses residual correlations (default: TRUE).
#'   When FALSE, the "rescor" prior class is automatically excluded even if listed in `what`.
#' @param use_empirical_legendre_priors Logical, use empirical priors for Legendre
#'   coefficients P1-P4 (default: TRUE). Applies timepoint-specific priors.
#' @param b_prior_sd Numeric, SD for fixed effect coefficients (default: NULL)
#' @param glucose_cv_residual Numeric, glucose residual CV = CVa + CVpre (default: 0.029)
#' @param glucose_cv_sd Numeric, SD for glucose residual CV prior (default: 1.5)
#' @param insulin_cv_residual Numeric, insulin residual CV = CVa + CVpre (default: 0.079)
#' @param insulin_cv_sd Numeric, scale for insulin residual CV prior (default: 1.5)
#' @param insulin_prior_df Numeric, degrees of freedom for Student-t insulin prior (default: 3)
#' @param sigma_slope_sd Numeric, SD for sigma slopes (default: 0.1)
#' @param rescor_lkj Numeric, LKJ prior for (2N)×(2N) residual correlation (default: 1.1)
#' @param L_lkj Numeric, LKJ prior for random effect correlation (default: 1.1)
#' @param prior_sd_scale Positive scalar applied as a multiplicative widening
#'   factor to every SD prior the function emits (default: 1.0 = no change).
#'   Specifically multiplies: the empirical Legendre P1-P4 SDs (and the
#'   `b_prior_sd` fallback), the P0 intercept SD scale `10 * sd(y)`,
#'   `glucose_cv_sd`, `insulin_cv_sd`, and the per-response RE-SD scale
#'   `sd(y)`. Means/locations, `sigma_slope_sd`, and the LKJ shapes
#'   (`rescor_lkj`, `L_lkj`) are intentionally NOT scaled. Use values > 1 for
#'   diffuse sensitivity refits (e.g., `prior_sd_scale = 5`).
#' @param data Data frame containing the response columns
#'   (`logglucose<suffix>`, `loginsulin<suffix>`). Required whenever any
#'   data-scaled prior is emitted — i.e., when `"sd"` is in `what` (RE-SD
#'   priors) or when `"beta"` is in `what` and `use_empirical_legendre_priors`
#'   is TRUE (P0 intercept priors). Used to compute each response's empirical
#'   SD for:
#'   \itemize{
#'     \item P0 intercepts (`leg_poly_0_ctrl`, `leg_poly_0_treat`):
#'           `normal(0, 10 * prior_sd_scale * sd(y))`
#'     \item Random-effect SDs: `student_t(4, 0, sd(y) * prior_sd_scale)`
#'           (half-t since SDs are positive)
#'   }
#'   Pass the same data frame you'll feed to `brms::brm()`.
#' @param independence Logical, indicates the formula was built with
#'   `independence = TRUE` (default: FALSE). When TRUE, the `"etacor"` class is
#'   automatically removed from `what` (with a message) because the independence
#'   formula contains no `class="L"` random-effect correlation matrix. The
#'   `class="sd"` priors are unchanged and still apply to each independent SD
#'   parameter brms emits.
#'
#' @return brmsprior object
#'
#' @details
#' Residual variance priors are centered on CVa + CVpre (known laboratory error),
#' with spread wide enough to accommodate CVi (biological variation) if data support it.
#' This conservative approach lets data inform whether single-OGTT uncertainty
#' exceeds analytical error alone.
#'
#' Prior location (CVa + CVpre):
#'   - Glucose: ~2.9% (CVa ~2.5% + CVpre ~1.5%)
#'   - Insulin: ~7.9% (CVa ~7.5% + CVpre ~2.5%)
#'
#' Prior spread (wide enough for CVi):
#'   - Glucose: Normal with SD=0.5; 95th percentile ~8% CV (covers CVi ~5.5%)
#'   - Insulin: Student-t(3) with scale=0.5; 95th percentile ~40% CV (covers CVi ~27%)
#'
#' Same analyte shares residual CV prior across all timepoints.
#' (2N)×(2N) rescor enables cross-timepoint and cross-analyte correlations.
#'
#' Prior parameter classes:
#'   - "beta": Fixed effects (Legendre polynomial coefficients)
#'   - "sigma": Residual standard deviation intercepts
#'   - "wiggle": Residual SD slopes (heteroskedasticity across timepoints)
#'   - "rescor": Residual correlations (between outcomes)
#'   - "sd": Random effect standard deviations
#'   - "etacor": Random effect correlations (latent variable correlations, L matrix)
#'
#' @examples
#' # 2 pregnancy timepoints (baseline 31wk + post 36wk) with all priors
#' priors_01 <- create_joined_growth_curve_priors_n(c("0", "1"), data = mydat)
#'
#' # Fixed effects model: no random effect priors
#' priors_fixed <- create_joined_growth_curve_priors_n(
#'   c("0", "1"),
#'   what = c("beta", "sigma", "wiggle", "rescor"),
#'   data = mydat
#' )
#'
#' # Prior sensitivity: beta priors only
#' priors_beta <- create_joined_growth_curve_priors_n(
#'   c("0", "1"), what = "beta", data = mydat
#' )
#'
#' # Prior sensitivity: weakly informative Legendre priors
#' priors_weak <- create_joined_growth_curve_priors_n(
#'   c("0", "1"),
#'   use_empirical_legendre_priors = FALSE,
#'   b_prior_sd = 5,
#'   data = mydat
#' )
#'
#' @export
create_joined_growth_curve_priors_n <- function(
    suffixes,
    what = c("beta", "sigma", "wiggle", "rescor", "sd", "etacor"),
    rescor = TRUE,
    use_empirical_legendre_priors = TRUE,
    b_prior_sd = NULL,
    glucose_cv_residual = 0.029,
    glucose_cv_sd = 1.5,
    insulin_cv_residual = 0.079,
    insulin_cv_sd = 1.5,
    insulin_prior_df = 3,
    sigma_slope_sd = 0.1,
    rescor_lkj = 1.1,
    L_lkj = 1.1,
    add_std_dev_lb = FALSE,
    prior_sd_scale = 1.0,
    data = NULL,
    independence = FALSE
) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' required but not installed.")
  }

  if (length(suffixes) < 2) {
    stop("Need at least 2 timepoints. Use single-timepoint priors for N=1.")
  }

  if (!is.numeric(prior_sd_scale) || length(prior_sd_scale) != 1L ||
      !is.finite(prior_sd_scale) || prior_sd_scale <= 0) {
    stop("prior_sd_scale must be a positive finite scalar.")
  }

  if (!is.logical(independence) || length(independence) != 1L || is.na(independence)) {
    stop("independence must be a single logical (TRUE/FALSE).")
  }

  # Validate what argument
  valid_what <- c("beta", "sigma", "wiggle", "rescor", "sd", "etacor")
  invalid <- setdiff(what, valid_what)
  if (length(invalid) > 0) {
    stop("Invalid 'what' values: ", paste(invalid, collapse = ", "),
         "\nValid options: ", paste(valid_what, collapse = ", "))
  }

  # Independence formula has no class="L" correlation matrix — drop etacor.
  if (independence && "etacor" %in% what) {
    message("independence = TRUE: removing 'etacor' from `what` (no class='L' correlation matrix to prior).")
    what <- setdiff(what, "etacor")
  }
  
  # Generate outcome variable names
  glucose_vars <- paste0("logglucose", suffixes)
  insulin_vars <- paste0("loginsulin", suffixes)
  all_outcomes <- c(glucose_vars, insulin_vars)

  # Per-response empirical SDs — required when "sd" is in `what` (RE-SD priors)
  # or when use_empirical_legendre_priors is TRUE and "beta" is in `what`
  # (P0 intercept priors are data-scaled).
  needs_data <- ("sd" %in% what) ||
                ("beta" %in% what && use_empirical_legendre_priors)
  emp_sd <- NULL
  if (needs_data) {
    if (is.null(data)) {
      stop("`data` must be provided: per-response data-scaled priors require ",
           "the same data frame used in brm() to compute sd(y) per outcome.")
    }
    missing_outcomes <- setdiff(all_outcomes, names(data))
    if (length(missing_outcomes) > 0) {
      stop("Missing outcome columns in `data`: ",
           paste(missing_outcomes, collapse = ", "))
    }
    emp_sd <- vapply(all_outcomes,
                     function(y) stats::sd(data[[y]], na.rm = TRUE),
                     numeric(1))
    bad_sd <- !is.finite(emp_sd) | emp_sd <= 0
    if (any(bad_sd)) {
      stop("Non-positive or non-finite empirical SD for: ",
           paste(all_outcomes[bad_sd], collapse = ", "))
    }
  }

  # Collect priors
  prior_list <- list()

  # Fixed effect coefficients (Legendre polynomial coefficients)
  if ("beta" %in% what) {
    if (use_empirical_legendre_priors) {
      # Use timepoint-specific empirical priors for P1-P4 for each suffix
      for (suffix in suffixes) {
        prior_list <- c(
          prior_list,
          list(create_legendre_coefficient_priors(suffix, "glucose", scale = prior_sd_scale)),
          list(create_legendre_coefficient_priors(suffix, "insulin", scale = prior_sd_scale))
        )
      }
      # P0 (leg_poly_0_ctrl, leg_poly_0_treat) intercept priors per response.
      # Data-scaled: normal(0, 10 * prior_sd_scale * sd(y)).
      for (outcome in all_outcomes) {
        p0_scale <- 10 * prior_sd_scale * emp_sd[[outcome]]
        for (grp in c("ctrl", "treat")) {
          prior_list <- c(prior_list, list(brms::set_prior(
            sprintf("normal(0, %s)", signif(p0_scale, 4)),
            class = "b",
            coef  = sprintf("leg_poly_0_%s", grp),
            resp  = outcome
          )))
        }
      }
    } else if (!is.null(b_prior_sd)) {
      # Fallback: uniform weakly informative prior for all coefficients
      b_prior_sd_scaled <- b_prior_sd * prior_sd_scale
      for (outcome in all_outcomes) {
        prior_list <- c(prior_list, list(
          brms::set_prior(
            paste0("normal(0, ", signif(b_prior_sd_scaled, 4), ")"),
            class = "b",
            resp = outcome
          )
        ))
      }
    }
  }
  
  # Glucose residual SD intercept (all glucose outcomes)
  # Note: This is on log scale (log link), so no lower bound needed
  if ("sigma" %in% what && !is.null(glucose_cv_residual)) {
    glucose_sd <- if (!is.null(glucose_cv_sd)) glucose_cv_sd else 1.5
    glucose_sd <- glucose_sd * prior_sd_scale
    for (outcome in glucose_vars) {
      prior_list <- c(prior_list, list(
        brms::set_prior(
          paste0("normal(log(sqrt(log(1 + ", glucose_cv_residual, "^2))), ",
                 signif(glucose_sd, 4), ")"),
          class = "Intercept",
          resp = outcome,
          dpar = "sigma"
        )
      ))
    }
  }
  
  # Insulin residual SD intercept (all insulin outcomes)
  # Uses Student-t(df=3) for fat tails to accommodate empirical extremes
  # Note: This is on log scale (log link), so no lower bound needed
  if ("sigma" %in% what && !is.null(insulin_cv_residual)) {
    insulin_scale <- if (!is.null(insulin_cv_sd)) insulin_cv_sd else 1.5
    insulin_scale <- insulin_scale * prior_sd_scale
    insulin_location <- paste0("log(sqrt(log(1 + ", insulin_cv_residual, "^2)))")

    # Use Student-t if df specified, otherwise Normal
    if (!is.null(insulin_prior_df)) {
      prior_str <- paste0("student_t(", insulin_prior_df, ", ", insulin_location, ", ", signif(insulin_scale, 4), ")")
    } else {
      prior_str <- paste0("normal(", insulin_location, ", ", signif(insulin_scale, 4), ")")
    }
    
    for (outcome in insulin_vars) {
      prior_list <- c(prior_list, list(
        brms::set_prior(
          prior_str,
          class = "Intercept",
          resp = outcome,
          dpar = "sigma"
        )
      ))
    }
  }
  
  # Residual SD slopes (all outcomes) - heteroskedasticity
  # Note: This is on log scale (log link), so no lower bound needed
  if ("wiggle" %in% what && !is.null(sigma_slope_sd)) {
    for (outcome in all_outcomes) {
      prior_list <- c(prior_list, list(
        brms::set_prior(
          paste0("normal(0, ", sigma_slope_sd, ")"),
          class = "b",
          resp = outcome,
          dpar = "sigma"
        )
      ))
    }
  }
  
  # Residual correlation ((2N)×(2N) matrix) — skip when rescor = FALSE
  if ("rescor" %in% what && !is.null(rescor_lkj) && rescor) {
    prior_list <- c(prior_list, list(
      brms::set_prior(
        paste0("lkj(", rescor_lkj, ")"),
        class = "rescor"
      )
    ))
  }
  
  # Cholesky factor for random effect correlation (eta correlations)
  if ("etacor" %in% what && !is.null(L_lkj)) {
    prior_list <- c(prior_list, list(
      brms::set_prior(
        paste0("lkj(", L_lkj, ")"),
        class = "L"
      )
    ))
  }
  
  # Random effect standard deviations (all outcomes)
  # Hardcoded data-scaled prior per response:
  #   student_t(4, 0, sd(y) * prior_sd_scale)
  # brms applies the implicit half-t truncation at 0 since SDs are positive.
  # `emp_sd` was validated/computed at the top of the function.
  if ("sd" %in% what) {
    for (outcome in all_outcomes) {
      sd_value <- emp_sd[[outcome]] * prior_sd_scale
      sd_prior_args <- list(
        prior = paste0("student_t(4, 0, ", signif(sd_value, 4), ")"),
        class = "sd",
        resp  = outcome
      )
      if (add_std_dev_lb) sd_prior_args$lb <- 1e-3
      prior_list <- c(prior_list, list(do.call(brms::set_prior, sd_prior_args)))
    }
  }
  
  # Combine all priors
  if (length(prior_list) == 0) {
    return(NULL)
  } else {
    return(do.call(c, prior_list))
  }
}

