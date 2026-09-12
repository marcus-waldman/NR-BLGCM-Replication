#' Create Formula for N-Timepoint Joined Growth Curve Model
#'
#' Generalized version of create_joined_growth_curve_formula() that handles
#' arbitrary number of timepoints. Generates formula for 2N outcomes
#' (N timepoints × 2 analytes).
#'
#' @param suffixes Character vector of timepoint suffixes (e.g., c("0", "1", "2"))
#' @param shared_random_effects Logical, share random effects across outcomes (default: TRUE)
#' @param effects Character, "random" (default) or "fixed". Controls individual-level structure:
#'   - "random": Latent variable model with random effects, multiple-group structure (default)
#'   - "fixed": Fixed effects for PTID with Legendre interactions, no grouping structure
#' @param rescor Logical, estimate residual correlations between outcomes (default: TRUE).
#'   Set to FALSE when glucose-insulin covariance is captured by shared random effects (| a |).
#' @param independence Logical, fit an independence random-effects model (default: FALSE).
#'   When TRUE, each outcome gets its own random-effects term using `||` syntax, with no
#'   shared `| a |` identifier across outcomes. This removes all correlations among random
#'   effects: cross-outcome, cross-timepoint, and within-outcome (P0..P4 are independent).
#'   Requires `effects = "random"`. Combine with `rescor = FALSE` for a fully independent
#'   baseline (no parameter anywhere models cross-outcome dependence).
#'
#' @return mvbrmsformula object with 2N outcomes
#'
#' @details
#' Random effects model (default):
#'   - 2N responses share identical Legendre predictors
#'   - Random effects: (leg_poly_0..4 | a | gr(PTID, by=treat)) shared across outcomes
#'   - Multiple-group structure with group-specific fixed effects (leg_poly_X_ctrl, leg_poly_X_treat)
#'   - Heteroskedastic: sigma ~ 1 + model_timepoint
#'   - rescor = TRUE enables (2N)×(2N) residual correlation
#'
#' Independence random effects model (independence = TRUE):
#'   - Per-outcome RE term: (leg_poly_0..4 || gr(PTID, by=treat)), no shared `a`
#'   - All 2N x 5 random effects are mutually uncorrelated
#'   - No class="L" correlation matrix is estimated
#'   - Drop the corresponding LKJ prior in the prior specification
#'
#' Fixed effects model:
#'   - Saturated individual-level model (no pooling)
#'   - Formula: y ~ 0 + PTID*leg_poly_0 + PTID*leg_poly_1 + ... + PTID*leg_poly_4
#'   - Each person gets explicit polynomial coefficients (no reference category)
#'   - NO multiple-group structure (would be redundant/non-identified with PTID fixed effects)
#'   - Useful for simulation validation (compare estimated vs. true individual parameters)
#'
#' Outcome ordering: All glucose outcomes, then all insulin outcomes
#'   (logglucose0, logglucose1, ..., loginsulin0, loginsulin1, ...)
#'
#' @examples
#' # 3 timepoints with random effects (default)
#' mvform_012 <- create_joined_growth_curve_formula_n(c("0", "1", "2"))
#'
#' # Fixed effects for simulation validation
#' mvform_012_fixed <- create_joined_growth_curve_formula_n(c("0", "1", "2"), effects = "fixed")
#'
#' # Independence baseline (all random effects independent, no rescor)
#' mvform_01_indep <- create_joined_growth_curve_formula_n(
#'   c("0", "1"), rescor = FALSE, independence = TRUE
#' )
#'
#' @export
create_joined_growth_curve_formula_n <- function(suffixes, shared_random_effects = TRUE, effects = "random", rescor = TRUE, independence = FALSE) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' required but not installed.")
  }

  if (length(suffixes) < 2) {
    stop("Need at least 2 timepoints. Use single-timepoint model for N=1.")
  }

  # Validate effects argument
  if (!effects %in% c("random", "fixed")) {
    stop("effects must be 'random' or 'fixed'")
  }

  if (!is.logical(independence) || length(independence) != 1L || is.na(independence)) {
    stop("independence must be a single logical (TRUE/FALSE).")
  }
  if (independence && effects != "random") {
    stop("independence = TRUE requires effects = 'random' (no random effects to decorrelate when effects = 'fixed').")
  }

  # Generate outcome variable names
  glucose_vars <- paste0("logglucose", suffixes)
  insulin_vars <- paste0("loginsulin", suffixes)

  if (effects == "random") {
    # Random effects specification.
    #  - Default ( | a | ): shared multivariate normal across all outcomes (single class="L").
    #  - Independence ( || ): per-outcome term with no shared identifier and within-term
    #    correlations dropped, so all REs are mutually independent.
    if (independence) {
      re_spec <- "(0 + leg_poly_0 + leg_poly_1 + leg_poly_2 + leg_poly_3 + leg_poly_4 || gr(PTID, by = treat))"
    } else {
      re_spec <- "(0 + leg_poly_0 + leg_poly_1 + leg_poly_2 + leg_poly_3 + leg_poly_4 | a | gr(PTID, by = treat))"
    }

    # Fixed effects specification (same for all outcomes)
    fixed_spec <- paste0(
      "0 + ",
      "leg_poly_0_ctrl + leg_poly_1_ctrl + leg_poly_2_ctrl + leg_poly_3_ctrl + leg_poly_4_ctrl + ",
      "leg_poly_0_treat + leg_poly_1_treat + leg_poly_2_treat + leg_poly_3_treat + leg_poly_4_treat + ",
      re_spec
    )
  } else {
    # Fixed effects specification (no grouping, no random effects)
    # Each individual gets their own polynomial explicitly estimated (no reference category)
    fixed_spec <- paste0(
      "0 + ",
      "PTID*leg_poly_0 + PTID*leg_poly_1 + PTID*leg_poly_2 + PTID*leg_poly_3 + PTID*leg_poly_4"
    )
  }
  
  # Create formulas for all glucose outcomes
  bf_glucose_list <- lapply(glucose_vars, function(var) {
    brms::bf(
      formula = as.formula(paste0(var, " | mi() ~ ", fixed_spec)),
      sigma ~ 1 + model_timepoint
    )
  })
  
  # Create formulas for all insulin outcomes
  bf_insulin_list <- lapply(insulin_vars, function(var) {
    brms::bf(
      formula = as.formula(paste0(var, " | mi() ~ ", fixed_spec)),
      sigma ~ 1 + model_timepoint
    )
  })
  
  # Combine all formulas
  all_bf_list <- c(bf_glucose_list, bf_insulin_list)
  
  # Build multivariate formula
  mvform <- do.call(brms::mvbrmsformula, all_bf_list) + brms::set_rescor(rescor)
  
  return(mvform)
}