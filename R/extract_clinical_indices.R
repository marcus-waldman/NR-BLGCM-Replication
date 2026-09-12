#' AUC treatment effect (cross-adjusted ANCOVA) at 36-week post-intervention.
#'
#' Posterior CHOICE - LC/Conventional treatment effect (estimate and 95%
#' interval) for glucose AUC and insulin AUC. For each posterior draw, per-PTID
#' AUCs are computed for both visits (suffix "0" = baseline 31wk, "1" = post
#' 36wk) and two ANCOVA regressions are fit:
#'
#'   AUC_glucose_1 ~ treat + AUC_glucose_0 + AUC_insulin_0
#'   AUC_insulin_1 ~ treat + AUC_glucose_0 + AUC_insulin_0
#'
#' Both baseline AUCs adjust each post outcome (cross-analyte ANCOVA). AUC is
#' trapezoidal over c(0, 30, 60, 90, 120) via bayestestR::area_under_curve;
#' in observed mode the t=-15 and t=0 fasting readings are averaged on the raw
#' scale and used as y[1].
#'
#' Modes:
#'   true_score - posterior_epred per response (alpha + zeta + Lambda baked
#'                in), per-draw ANCOVA, posterior summary of the treat coef.
#'   observed   - raw + Ymi imputations per response, per-draw ANCOVA, Rubin's
#'                rules pooling via mice::pool on a mira of per-draw lm fits.
#'
#' Requires legendre_polynomials.R sourced (leg_poly, scale_to_legendre_domain).

# Trapezoidal AUC over the five clinical OGTT timepoints (Δt = 30). For x =
# c(0, 30, 60, 90, 120), the rule reduces to a linear combination of y with
# weights .AUC_TRAP_W = c(15, 30, 30, 30, 15) (0.5 * Δt * adjacency counts),
# so vector AUC is sum(w * y) and matrix AUC is a single Y %*% w. This is
# mathematically identical to bayestestR::area_under_curve on this grid and
# matches it to machine precision (verified in tmp/demo_auc_vectorize.R).
.AUC_TRAP_W <- c(15, 30, 30, 30, 15)

# Polymorphic: y can be a length-5 vector (returns scalar) or an (n x 5) matrix
# of n curves (returns length-n vector of AUCs). Existing scalar callers see
# unchanged behaviour; new vectorised call sites pass a matrix.
auc_trap <- function(y) {
  if (is.matrix(y) || is.array(y) && length(dim(y)) == 2) {
    return(as.numeric(y %*% .AUC_TRAP_W))
  }
  return(sum(.AUC_TRAP_W * y))
}

# Construct the newdata that posterior_epred expects: every PTID crossed with
# the 5 clinical timepoints (46 x 5 = 230 rows), with the Legendre basis
# pre-computed and split into the _ctrl / _treat interaction columns referenced
# by the model formula (see create_joined_growth_curve_formula_n). The same
# newdata works for all four responses (logglucose0/1, loginsulin0/1) because
# the structural design is shared across them.
build_clinical_newdata <- function(fit) {
  #Verified
  fit$data |>
    dplyr::distinct(PTID, treat) |>                  # 46 x 2 (PTID, treat)
    dplyr::arrange(PTID) |>
    tidyr::crossing(t = c(0, 30, 60, 90, 120)) |>    # 46 x 5 = 230 rows
    dplyr::mutate(
      t_std            = scale_to_legendre_domain(t),
      leg_poly_0       = 1,
      leg_poly_1       = leg_poly(t_std, 1),
      leg_poly_2       = leg_poly(t_std, 2),
      leg_poly_3       = leg_poly(t_std, 3),
      leg_poly_4       = leg_poly(t_std, 4),
      leg_poly_0_ctrl  = leg_poly_0 * (1 - treat),
      leg_poly_1_ctrl  = leg_poly_1 * (1 - treat),
      leg_poly_2_ctrl  = leg_poly_2 * (1 - treat),
      leg_poly_3_ctrl  = leg_poly_3 * (1 - treat),
      leg_poly_4_ctrl  = leg_poly_4 * (1 - treat),
      leg_poly_0_treat = leg_poly_0 * treat,
      leg_poly_1_treat = leg_poly_1 * treat,
      leg_poly_2_treat = leg_poly_2 * treat,
      leg_poly_3_treat = leg_poly_3 * treat,
      leg_poly_4_treat = leg_poly_4 * treat,
      model_timepoint  = factor(t, levels = c("0", "30", "60", "90", "120"))
    ) |>
    dplyr::select(-t) |>
    dplyr::arrange(PTID, model_timepoint) %>%
    dplyr::relocate(model_timepoint, .after = "PTID")
}

# Per-draw, per-PTID AUC matrix for one response, in true-score mode.
# Returns [n_draws x n_indiv]. Uses the model's structural expectation
# (posterior_epred), so PTID random effects are integrated in per draw.
true_score_auc_one_resp <- function(fit, nd, resp, n_indiv) {
  ep <- brms::posterior_epred(fit, newdata = nd, resp = resp)
  n_draws <- nrow(ep)
  # nd is in (PTID, time) order; column-major reshape puts time on axis 2,
  # PTID on axis 3. Replace apply(..., c(1,3), auc_trap) -- the previous hot
  # path -- with the explicit five-term weighted sum on slabs of the cube. Each
  # slab ep5[, t, ] is an (n_draws x n_indiv) matrix, so this is five
  # element-wise scalar-times-matrix products and four matrix adds: fully
  # vectorised, no apply, no per-cell call into bayestestR.
  ep5 <- exp(ep)
  dim(ep5) <- c(n_draws, 5, n_indiv)
  w <- .AUC_TRAP_W
  return(w[1] * ep5[, 1, ] + w[2] * ep5[, 2, ] + w[3] * ep5[, 3, ] +
         w[4] * ep5[, 4, ] + w[5] * ep5[, 5, ])
}

# Fit per-draw cross-adjusted ANCOVA for one post outcome (glucose or insulin)
# and return the per-draw treat coefficient, optionally with its standard
# error for downstream Rubin's-rules pooling.
#
# Dispatch by formula:
#   * Unadjusted (RHS = "treat"): X is fixed across draws. A single
#     multi-response .lm.fit(X, Y_matrix) call replaces n_draws per-draw
#     calls (~9x faster); SE uses the fit's own QR so (tau, se) is bytewise
#     identical to summary(lm(...))$coefficients.
#   * Adjusted (RHS = treat + AUC_glucose_0 + AUC_insulin_0): X varies per
#     draw, so we loop with .lm.fit + chol2inv(qr.R(fit$qr)) for the SE.
#     Bytewise identical to summary(lm(...))$coefficients.
# Returns:
#   return_fits = FALSE  -> numeric(n_draws) of per-draw tau
#   return_fits = TRUE   -> list(tau = num, se = num, dfcom = scalar)
ancova_treat_per_draw <- function(formula, auc_glu_0, auc_ins_0, auc_glu_1, auc_ins_1,
                                  treat_vec, return_fits = FALSE) {
  rhs <- attr(stats::terms(formula), "term.labels")
  lhs <- all.vars(formula)[1]
  y_mat <- switch(lhs,
    "AUC_glucose_1" = auc_glu_1,
    "AUC_insulin_1" = auc_ins_1,
    stop("Unsupported LHS: ", lhs)
  )
  if (identical(rhs, "treat")) {
    return(.ancova_unadj(y_mat, treat_vec, return_fits))
  }
  if (setequal(rhs, c("treat", "AUC_glucose_0", "AUC_insulin_0"))) {
    return(.ancova_adj(y_mat, auc_glu_0, auc_ins_0, treat_vec, return_fits))
  }
  stop("Unsupported RHS: ", paste(rhs, collapse = " + "))
}

# Unadjusted ANCOVA: X = [1, treat] fixed across draws, so a single
# multi-response .lm.fit handles all n_draws regressions at once.
.ancova_unadj <- function(y_mat, treat_vec, return_fits) {
  n_indiv <- length(treat_vec)
  X <- cbind(1, treat_vec)
  Y <- t(y_mat)                                       # n_indiv x n_draws
  fit <- .lm.fit(X, Y)
  tau <- as.numeric(fit$coefficients[2, ])
  if (!return_fits) { return(tau) }
  dfcom <- n_indiv - 2
  sigma2 <- colSums(fit$residuals^2) / dfcom          # length n_draws
  qr_obj <- structure(list(qr = fit$qr, rank = fit$rank,
                           qraux = fit$qraux, pivot = fit$pivot),
                      class = "qr")
  xtx_inv_22 <- chol2inv(qr.R(qr_obj))[2, 2]
  se <- sqrt(sigma2 * xtx_inv_22)
  return(list(tau = tau, se = se, dfcom = dfcom))
}

# Adjusted ANCOVA: X = [1, treat, AUC_glu_0[d,], AUC_ins_0[d,]] varies per
# draw. Stays on .lm.fit per draw because (a) a batched-normal-equations
# variant we tried (vectorised X'X cross-products + per-draw 4x4 solve) was
# no faster than .lm.fit per draw at this scale -- the loop overhead
# dominates either way -- and (b) directly forming X'X and solving it
# squares the condition number, introducing ~1e-9 drift in tau on our
# AUC-scale data; QR-based OLS via .lm.fit is bytewise to summary.lm.
.ancova_adj <- function(y_mat, auc_glu_0, auc_ins_0, treat_vec, return_fits) {
  n_indiv <- length(treat_vec)
  n_draws <- nrow(y_mat)
  dfcom <- n_indiv - 4
  intercept <- rep(1, n_indiv)
  tau <- numeric(n_draws)
  if (!return_fits) {
    for (d in seq_len(n_draws)) {
      X <- cbind(intercept, treat_vec, auc_glu_0[d, ], auc_ins_0[d, ])
      tau[d] <- .lm.fit(X, y_mat[d, ])$coefficients[2]
    }
    return(tau)
  }
  se <- numeric(n_draws)
  for (d in seq_len(n_draws)) {
    X <- cbind(intercept, treat_vec, auc_glu_0[d, ], auc_ins_0[d, ])
    fit <- .lm.fit(X, y_mat[d, ])
    tau[d] <- fit$coefficients[2]
    sigma2 <- sum(fit$residuals^2) / dfcom
    qr_obj <- structure(list(qr = fit$qr, rank = fit$rank,
                             qraux = fit$qraux, pivot = fit$pivot),
                        class = "qr")
    se[d] <- sqrt(sigma2 * chol2inv(qr.R(qr_obj))[2, 2])
  }
  return(list(tau = tau, se = se, dfcom = dfcom))
}

# True-score mode: posterior_epred returns Lambda %*% (alpha + zeta_i) per
# draw with PTID random effects baked in. Per draw, fit cross-adjusted ANCOVA
# on the four per-PTID AUCs and summarize the treat coefficient as posterior
# mean + 95% CrI. No within-draw SE; no Rubin's pooling.
compute_auc_true_score <- function(fit) {
  nd <- build_clinical_newdata(fit)
  ptid_info <- dplyr::distinct(nd, PTID, treat) |> dplyr::arrange(PTID)
  n_indiv <- nrow(ptid_info)

  auc_glu_0 <- true_score_auc_one_resp(fit, nd, "logglucose0", n_indiv)
  auc_ins_0 <- true_score_auc_one_resp(fit, nd, "loginsulin0", n_indiv)
  auc_glu_1 <- true_score_auc_one_resp(fit, nd, "logglucose1", n_indiv)
  auc_ins_1 <- true_score_auc_one_resp(fit, nd, "loginsulin1", n_indiv)

  tau_glu <- ancova_treat_per_draw(
    AUC_glucose_1 ~ treat + AUC_glucose_0 + AUC_insulin_0,
    auc_glu_0, auc_ins_0, auc_glu_1, auc_ins_1, ptid_info$treat
  )
  tau_ins <- ancova_treat_per_draw(
    AUC_insulin_1 ~ treat + AUC_glucose_0 + AUC_insulin_0,
    auc_glu_0, auc_ins_0, auc_glu_1, auc_ins_1, ptid_info$treat
  )

  list(
    AUC_glucose = list(estimate = mean(tau_glu), ci = unname(quantile(tau_glu, c(0.025, 0.975)))),
    AUC_insulin = list(estimate = mean(tau_ins), ci = unname(quantile(tau_ins, c(0.025, 0.975))))
  )
}

# brms uses one of two bracket-index conventions for Ymi columns: sequential
# 1..N_mi, or data-row indices (= jmi). We detect which is in play and return
# column names ordered so that the k-th entry corresponds to the k-th missing
# observation (whose row in fit$data is jmi[k]).
match_ymi_cols <- function(post, prefix, jmi) {
  #Validated
  cols <- grep(paste0("^", prefix, "\\["), names(post), value = TRUE)
  if (length(cols) == 0) stop("No '", prefix, "[k]' columns in posterior")

  brk <- as.integer(sub(paste0("^", prefix, "\\[(\\d+)\\]$"), "\\1", cols))
  lookup <- setNames(cols, brk)

  if (setequal(brk, jmi))             lookup[as.character(jmi)]
  else if (setequal(brk, seq_along(jmi))) lookup[as.character(seq_along(jmi))]
  else stop(prefix, " bracket indices match neither Jmi nor 1..", length(jmi))
}

# Reshape model_data so that rows are ordered (PTID, time) with the two
# fasting rows (t=-15 then t=0) in that order within each PTID. The
# model_timepoint factor only has 5 levels (c("0","30","60","90","120"));
# both fasting rows carry the "0" level, so cumsum(model_timepoint=="0")
# breaks the tie (1 -> t=-15, 2 -> t=0). Tags .row = original row index in
# fit$data so we can later place Ymi[k] into obs_data via match(jmi, .row).
prepare_obs_data <- function(model_data, n_indiv) {
  obs_data <- model_data |>
    dplyr::mutate(.row = seq_len(dplyr::n())) |>
    dplyr::group_by(PTID) |>
    dplyr::mutate(.t0 = cumsum(model_timepoint == "0")) |>
    dplyr::ungroup() |>
    dplyr::arrange(PTID, model_timepoint, .t0) |>
    dplyr::select(-.t0)
  stopifnot(nrow(obs_data) == n_indiv * 6)
  obs_data
}

# Per-draw, per-PTID AUC matrix for one response, in observed mode. For each
# draw, missing log-values are filled from that draw's Ymi imputations, the
# two fasting columns are averaged on the raw scale, and trapezoidal AUC is
# computed at (0, 30, 60, 90, 120). Returns [n_draws x n_indiv]. Handles the
# zero-missingness case (Jmi empty -> no fill).
observed_auc_one_resp <- function(obs_data, post, sd, resp, n_indiv) {
  jmi <- sd[[paste0("Jmi_", resp)]]
  if (is.null(jmi)) jmi <- integer(0)

  if (length(jmi) > 0) {
    ymi_mat <- as.matrix(post[, match_ymi_cols(post, paste0("Ymi_", resp), jmi), drop = FALSE])
    fill_pos <- match(jmi, obs_data$.row)
    stopifnot(!anyNA(fill_pos))
  } else {
    ymi_mat <- NULL
    fill_pos <- integer(0)
  }

  n_draws <- nrow(post)
  base_log <- obs_data[[resp]]
  n_obs <- length(base_log)               # = 6 * n_indiv
  # Fully vectorised replacement for the original per-draw loop. Build the
  # (n_draws x n_obs) matrix of log-scale outcomes by broadcasting base_log
  # and overwriting the missing-position columns with the per-draw Ymi draws.
  # Reshape to (n_draws, 6, n_indiv) -- column-major puts the within-PTID time
  # axis on dim 2 -- exponentiate, then take the trapezoidal AUC with the t=1
  # cell averaged from the t=-15 and fasting readings (cols 1,2), per the
  # observed-mode definition. Final form is six element-wise scalar*matrix
  # products plus adds: no for loop, no apply.
  y_log_all <- matrix(base_log, nrow = n_draws, ncol = n_obs, byrow = TRUE)
  if (length(fill_pos)) {
    y_log_all[, fill_pos] <- ymi_mat
  }
  y6_all <- exp(y_log_all)
  dim(y6_all) <- c(n_draws, 6, n_indiv)
  w <- .AUC_TRAP_W
  y5_t1 <- 0.5 * (y6_all[, 1, ] + y6_all[, 2, ])
  return(w[1] * y5_t1 +
         w[2] * y6_all[, 3, ] + w[3] * y6_all[, 4, ] +
         w[4] * y6_all[, 5, ] + w[5] * y6_all[, 6, ])
}

# Pool per-draw lm fits via Rubin's rules. Wraps the fit list as a mira object
# so mice::pool() does within/between variance combination on the treat
# regression coefficient directly (which is the right within-imputation
# variance for an ANCOVA estimator, unlike the pooled two-sample SE used by
# the previous mean-difference version).
pool_treat_rubins <- function(res, conf_level = 0.95) {
  # Hand-rolled Rubin's-rules pooler for a single coefficient ("treat"),
  # bytewise equivalent to summary(mice::pool(mira_obj), conf.int = TRUE) on
  # the treat row. Replicates mice 3.18.0 mice:::pool.vector (rule
  # "rubin1987") + mice:::barnard.rubin, verified to 0.000e+00 against
  # mice::pool on synthetic per-draw lm fits.
  #
  # res: list(tau = numeric[m], se = numeric[m], dfcom = scalar) as returned
  # by ancova_treat_per_draw(..., return_fits = TRUE).
  est <- res$tau
  se  <- res$se
  m   <- length(est)
  qbar  <- mean(est)
  ubar  <- mean(se^2)
  b     <- stats::var(est)
  t_var <- ubar + (1 + 1/m) * b
  # mice:::barnard.rubin: clamp lambda below 1e-4 to avoid divide-by-zero in
  # dfold; pass dfcom = Inf through to skip the Barnard-Rubin observed-df arm.
  lambda <- (1 + 1/m) * b / t_var
  if (lambda < 1e-04) { lambda <- 1e-04 }
  dfold <- (m - 1) / lambda^2
  dfobs <- (res$dfcom + 1) / (res$dfcom + 3) * res$dfcom * (1 - lambda)
  df_BR <- if (is.infinite(res$dfcom)) { dfold } else { dfold * dfobs / (dfold + dfobs) }
  se_p  <- sqrt(t_var)
  tq    <- stats::qt((1 + conf_level) / 2, df_BR)
  return(list(estimate = qbar, ci = c(qbar - tq * se_p, qbar + tq * se_p)))
}

# Observed mode: for each posterior draw, fill the four responses' missing
# log-values from Ymi[draw], compute per-PTID AUCs at baseline (suffix 0) and
# post (suffix 1), fit cross-adjusted ANCOVA for each post outcome, and pool
# the treat coefficient across draws via Rubin's rules (each draw = one
# completed-data analysis).
compute_auc_observed <- function(fit) {
  model_data <- fit$data
  ptid_info  <- dplyr::distinct(model_data, PTID, treat) |> dplyr::arrange(PTID)
  n_indiv    <- nrow(ptid_info)

  obs_data <- prepare_obs_data(model_data, n_indiv)
  sd       <- brms::standata(fit)
  post     <- as.data.frame(brms::as_draws_df(fit))

  auc_glu_0 <- observed_auc_one_resp(obs_data, post, sd, "logglucose0", n_indiv)
  auc_ins_0 <- observed_auc_one_resp(obs_data, post, sd, "loginsulin0", n_indiv)
  auc_glu_1 <- observed_auc_one_resp(obs_data, post, sd, "logglucose1", n_indiv)
  auc_ins_1 <- observed_auc_one_resp(obs_data, post, sd, "loginsulin1", n_indiv)

  res_glu <- ancova_treat_per_draw(
    AUC_glucose_1 ~ treat + AUC_glucose_0 + AUC_insulin_0,
    auc_glu_0, auc_ins_0, auc_glu_1, auc_ins_1, ptid_info$treat,
    return_fits = TRUE
  )
  res_ins <- ancova_treat_per_draw(
    AUC_insulin_1 ~ treat + AUC_glucose_0 + AUC_insulin_0,
    auc_glu_0, auc_ins_0, auc_glu_1, auc_ins_1, ptid_info$treat,
    return_fits = TRUE
  )

  list(
    AUC_glucose = pool_treat_rubins(res_glu),
    AUC_insulin = pool_treat_rubins(res_ins)
  )
}

#' @param fit  brmsfit (joined_01).
#' @param mode "true_score" or "observed".
#' @return list(AUC_glucose, AUC_insulin), each list(estimate, ci).
#' @export
extract_auc_treatment_effect <- function(fit, mode = c("true_score", "observed")) {
  mode <- match.arg(mode)
  if (!inherits(fit, "brmsfit")) stop("fit must be a brmsfit object")
  if (!exists("leg_poly", mode = "function") ||
      !exists("scale_to_legendre_domain", mode = "function")) {
    stop("Source R/legendre_polynomials.R first.")
  }
  for (pkg in c("brms", "bayestestR", "mice", "tidyr", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.")
  }
  if (mode == "true_score") compute_auc_true_score(fit) else compute_auc_observed(fit)
}
