#' OGTT Trajectory Summary Computation
#'
#' Functions for computing group-level trajectory summaries from brms posterior
#' draws. Supports both true-score smooth curves (individual-level marginalization
#' on a fine time grid) and observed group means (Rubin's combining rules at
#' clinical timepoints with Ymi imputation filling).
#'
#' @details
#' Both functions return data.frames with identical column structure:
#'   minute, mean, lo, hi, group, analyte
#' making them directly compatible for overlay plotting.
#'
#' Depends on: legendre_polynomials.R (legendre_design_matrix)


# ============================================================================
# TRUE-SCORE SMOOTH TRAJECTORIES
# ============================================================================

#' Compute True-Score Group Trajectories on Fine Time Grid
#'
#' For each posterior draw, computes individual-level trajectories
#' (eta_i = alpha_group + zeta_i) on a fine time grid, exponentiates to the
#' raw scale, averages within group, then summarizes across draws.
#'
#' @param fit A brms model fit object
#' @param suffixes Character vector of suffixes to include (e.g., "1" for post only)
#' @param analyte_specs Named list of analyte specifications. Each element is a list
#'   with `analyte` (display name) and `outcome` (brms outcome column). Default NULL
#'   auto-infers from suffix: logglucoseX, loginsulinX.
#' @param fine_times Numeric vector of timepoints for the fine grid (default seq(0, 120, by = 1))
#' @param max_degree Integer, maximum Legendre polynomial degree (default 4)
#' @param group_labels Named character vector mapping treat values to labels
#'   (default c("0" = "LC/Conventional", "1" = "High Carb"))
#' @param prob Numeric probability for credible intervals (default 0.95)
#'
#' @return data.frame with columns: minute, mean, lo, hi, group, analyte
#'
#' @details
#' The individual-level marginalization approach:
#'   1. Extract fixed-effect alpha (per group) and random-effect zeta (per individual)
#'   2. For each draw: eta_i = alpha_group + zeta_i -> Lambda_fine %*% eta_i -> exp()
#'   3. Average within group across individuals for each draw
#'   4. Summarize posterior mean and CrI across draws
#'
#' This produces smooth trajectories at 121 evaluation points (0 to 120 min).
#'
#' @export
compute_true_score_trajectories <- function(fit, suffixes = "1",
                                            analyte_specs = NULL,
                                            fine_times = seq(0, 120, by = 1),
                                            max_degree = 4,
                                            group_labels = c("0" = "LC/Conventional",
                                                             "1" = "CHOICE"),
                                            prob = 0.95) {
  if (!inherits(fit, "brmsfit")) {
    stop("fit must be a brmsfit object")
  }
  if (!exists("legendre_design_matrix", mode = "function")) {
    stop("Function legendre_design_matrix() not found. Source 'R/legendre_polynomials.R' first.")
  }
  
  alpha <- (1 - prob) / 2
  
  # Auto-infer analyte specs from suffix
  if (is.null(analyte_specs)) {
    analyte_specs <- lapply(suffixes, function(suf) {
      list(
        list(analyte = "Glucose", outcome = paste0("logglucose", suf)),
        list(analyte = "Insulin", outcome = paste0("loginsulin", suf))
      )
    })
    names(analyte_specs) <- suffixes
  } else {
    # If single list of specs provided, replicate for all suffixes
    if (!is.null(analyte_specs[[1]]$analyte)) {
      single_specs <- analyte_specs
      analyte_specs <- setNames(
        replicate(length(suffixes), single_specs, simplify = FALSE),
        suffixes
      )
    }
  }
  
  # Legendre design matrix on fine grid
  Lambda_fine <- legendre_design_matrix(fine_times, max_degree = max_degree)
  
  # Posterior draws
  
  draws <- brms::as_draws_df(fit)
  n_posterior <- nrow(draws)
  
  # Individual info
  model_data <- fit$data
  ptid_info <- model_data |>
    dplyr::select(PTID, treat) |>
    dplyr::distinct() |>
    dplyr::arrange(PTID)
  ptids <- as.character(ptid_info$PTID)
  n_individuals <- length(ptids)
  
  traj_data <- data.frame()
  
  for (suf in suffixes) {
    specs <- analyte_specs[[suf]]
    
    for (spec in specs) {
      cat("  ", spec$analyte, "(suffix ", suf, ")...\n", sep = "")
      
      # Pre-extract fixed effect draws for both groups: [n_posterior x (max_degree+1)]
      alpha_ctrl <- as.matrix(
        draws[, paste0("b_", spec$outcome, "_leg_poly_", 0:max_degree, "_ctrl")]
      )
      alpha_treat <- as.matrix(
        draws[, paste0("b_", spec$outcome, "_leg_poly_", 0:max_degree, "_treat")]
      )
      
      # Pre-extract random effects for all individuals
      ranef_list <- vector("list", n_individuals)
      for (i in seq_along(ptids)) {
        re_cols <- paste0("r_PTID__", spec$outcome, "[", ptids[i], ",leg_poly_", 0:max_degree, "]")
        ranef_list[[i]] <- as.matrix(draws[, re_cols])
      }
      
      # For each group, compute per-draw mean trajectory on raw scale
      for (grp_val_str in names(group_labels)) {
        grp_val <- as.numeric(grp_val_str)
        grp_label <- group_labels[grp_val_str]
        grp_alpha <- if (grp_val == 0) alpha_ctrl else alpha_treat
        grp_idx <- which(ptid_info$treat == grp_val)
        n_grp <- length(grp_idx)
        
        if (n_grp == 0) next
        
        # group_mean_traj[d, t]: mean raw-scale trajectory at timepoint t for draw d
        group_mean_traj <- matrix(0, nrow = n_posterior, ncol = length(fine_times))
        
        for (i in grp_idx) {
          eta_i <- grp_alpha + ranef_list[[i]]
          log_traj_i <- eta_i %*% t(Lambda_fine)
          group_mean_traj <- group_mean_traj + exp(log_traj_i)
        }
        group_mean_traj <- group_mean_traj / n_grp
        
        # Summarize across draws
        traj_mean <- colMeans(group_mean_traj)
        traj_lo <- apply(group_mean_traj, 2, quantile, probs = alpha)
        traj_hi <- apply(group_mean_traj, 2, quantile, probs = 1 - alpha)
        
        traj_data <- rbind(traj_data, data.frame(
          minute = fine_times,
          mean = traj_mean,
          lo = traj_lo,
          hi = traj_hi,
          group = grp_label,
          analyte = spec$analyte,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  rownames(traj_data) <- NULL
  traj_data
}


# ============================================================================
# OBSERVED GROUP MEANS (RUBIN'S COMBINING RULES)
# ============================================================================

#' Compute Observed Group Means at Clinical Timepoints
#'
#' For each group × analyte × clinical timepoint, fills missing values from
#' Ymi posterior imputations, computes per-draw cell means on the raw scale,
#' then combines using Rubin's rules for proper MI inference.
#'
#' @param fit A brms model fit object
#' @param suffixes Character vector of suffixes to include (e.g., "1" for post only)
#' @param analyte_specs Named list of analyte specifications. Each element is a list
#'   with `analyte` (display name) and `col` (brms column name). Default NULL
#'   auto-infers from suffix.
#' @param clinical_minutes Numeric vector of clinical timepoints (default c(0, 30, 60, 90, 120))
#' @param group_labels Named character vector mapping treat values to labels
#' @param n_within_subsample Integer, max draws for within-imputation variance estimation
#'   (default 200; subsampling for speed)
#'
#' @return data.frame with columns: minute, mean, lo, hi, group, analyte
#'
#' @details
#' For each cell (group × timepoint):
#'   1. Get raw observed values; fill NAs from Ymi[draw]
#'   2. Compute cell mean (exp scale) for each draw
#'   3. Pool using Rubin's combining rules:
#'      - Q_bar = mean of draw means
#'      - B = between-imputation variance
#'      - U_bar = average within-imputation variance
#'      - Total variance = U_bar + (1 + 1/m) * B
#'   4. Barnard-Rubin adjusted df for t-intervals
#'
#' @export
compute_observed_group_means <- function(fit, suffixes = "1",
                                         analyte_specs = NULL,
                                         clinical_minutes = c(0, 30, 60, 90, 120),
                                         group_labels = c("0" = "LC/Conventional",
                                                          "1" = "CHOICE"),
                                         n_within_subsample = 200) {
  if (!inherits(fit, "brmsfit")) {
    stop("fit must be a brmsfit object")
  }
  
  draws <- brms::as_draws_df(fit)
  n_posterior <- nrow(draws)
  model_data <- fit$data
  
  # Auto-infer analyte specs
  if (is.null(analyte_specs)) {
    analyte_specs <- lapply(suffixes, function(suf) {
      list(
        list(analyte = "Glucose", col = paste0("logglucose", suf)),
        list(analyte = "Insulin", col = paste0("loginsulin", suf))
      )
    })
    names(analyte_specs) <- suffixes
  } else {
    if (!is.null(analyte_specs[[1]]$analyte)) {
      single_specs <- analyte_specs
      analyte_specs <- setNames(
        replicate(length(suffixes), single_specs, simplify = FALSE),
        suffixes
      )
    }
  }
  
  obs_data <- data.frame()
  
  for (suf in suffixes) {
    specs <- analyte_specs[[suf]]
    
    for (spec in specs) {
      y_col <- spec$col
      
      # Get observed values and identify missing
      y_obs <- model_data[[y_col]]
      is_miss <- is.na(y_obs)
      
      # For each group × timepoint, compute per-draw means then pool
      for (grp_val_str in names(group_labels)) {
        grp_val <- as.numeric(grp_val_str)
        grp_label <- group_labels[grp_val_str]
        grp_mask <- model_data$treat == grp_val
        
        for (t_min in clinical_minutes) {
          t_mask <- model_data$model_timepoint == as.character(t_min)
          row_mask <- grp_mask & t_mask
          
          if (sum(row_mask) == 0) next
          
          # For each draw, fill Ymi and compute group mean at this timepoint
          y_base <- y_obs[row_mask]
          miss_in_subset <- is_miss[row_mask]
          global_indices <- which(row_mask)
          
          draw_means <- numeric(n_posterior)
          for (d in seq_len(n_posterior)) {
            y_filled <- y_base
            if (any(miss_in_subset)) {
              for (j in which(miss_in_subset)) {
                global_idx <- global_indices[j]
                ymi_col <- paste0("Ymi_", y_col, "[", global_idx, "]")
                if (ymi_col %in% names(draws)) {
                  y_filled[j] <- draws[[ymi_col]][d]
                }
              }
            }
            draw_means[d] <- mean(exp(y_filled), na.rm = TRUE)
          }
          
          # Rubin's combining rules
          Q_bar <- mean(draw_means)
          B <- var(draw_means)
          
          # Within-imputation variance (subsample for speed)
          n_in_cell <- sum(row_mask)
          n_sub <- min(n_posterior, n_within_subsample)
          U_bar <- 0
          for (d in seq_len(n_sub)) {
            y_filled <- y_base
            if (any(miss_in_subset)) {
              for (j in which(miss_in_subset)) {
                global_idx <- global_indices[j]
                ymi_col <- paste0("Ymi_", y_col, "[", global_idx, "]")
                if (ymi_col %in% names(draws)) {
                  y_filled[j] <- draws[[ymi_col]][d]
                }
              }
            }
            U_bar <- U_bar + var(exp(y_filled), na.rm = TRUE) / n_in_cell
          }
          U_bar <- U_bar / n_sub
          
          # Total variance (Rubin's rules)
          T_var <- U_bar + (1 + 1 / n_posterior) * B
          se <- sqrt(T_var)
          
          # Barnard-Rubin adjusted df
          m <- n_posterior
          df_num <- (m - 1) * (1 + U_bar / ((1 + 1/m) * B))^2
          df_use <- min(df_num, n_in_cell - 1)
          if (!is.finite(df_use) || df_use < 1) df_use <- n_in_cell - 1
          t_crit <- qt(0.975, df_use)
          
          obs_data <- rbind(obs_data, data.frame(
            minute = t_min,
            mean = Q_bar,
            lo = Q_bar - t_crit * se,
            hi = Q_bar + t_crit * se,
            group = grp_label,
            analyte = spec$analyte,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  rownames(obs_data) <- NULL
  obs_data
}
