# ============================================================================
# COMPOSED TRAJECTORY FIGURE
# ============================================================================

#' Compose Multi-Panel Trajectory Figure with Insets
#'
#' Assembles a two-panel (glucose + insulin) trajectory plot with AUC treatment
#' effect insets, using patchwork for layout.
#'
#' @param traj_data data.frame from compute_true_score_trajectories()
#' @param obs_data data.frame from compute_observed_group_means()
#' @param treat_results data.frame from compare_treatment_effects()
#' @param panel_specs Named list of panel specifications. Each element is a list
#'   with `analyte` (filter value), `y_label`, and `auc_label` (for inset filtering).
#'   Default provides Glucose and Insulin panels.
#' @param group_colours Named character vector of group colours
#' @param group_shapes Named integer vector of group point shapes
#' @param inset_position Numeric vector c(left, bottom, right, top) for inset placement
#'   (default c(0.35, 0.02, 0.98, 0.42))
#'
#' @return A patchwork object (composed ggplot)
#'
#' @export
compose_trajectory_figure <- function(traj_data, obs_data, treat_results,
                                      panel_specs = NULL,
                                      group_colours = c("LC/Conventional" = "#DF8F44",
                                                        "CHOICE" = "#00A1D5"),
                                      group_shapes = c("LC/Conventional" = 16,
                                                       "CHOICE" = 17),
                                      inset_position = c(0.35, 0.02, 0.98, 0.42),
                                      add_annotations = FALSE) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' required.")
  }
  
  # Default panel specs
  if (is.null(panel_specs)) {
    panel_specs <- list(
      glucose = list(analyte = "Glucose", y_label = "Glucose (mg/dL)",
                     auc_label = "Glucose AUC"),
      insulin = list(analyte = "Insulin", y_label = "Insulin (\u00B5IU/mL)",
                     auc_label = "Insulin AUC")
    )
  }
  
  panels <- list()
  for (spec_name in names(panel_specs)) {
    spec <- panel_specs[[spec_name]]
    
    traj_sub <- traj_data[traj_data$analyte == spec$analyte, ]
    obs_sub <- obs_data[obs_data$analyte == spec$analyte, ]
    
    p <- plot_trajectory_panel(traj_sub, obs_sub, spec$y_label,
                               group_colours, group_shapes)
    
    # Add annotations to first panel only (typically glucose)
    if (add_annotations && spec_name == names(panel_specs)[1]) {
      grp1 <- unique(traj_sub$group)[1]
      traj_grp1 <- traj_sub[traj_sub$group == grp1, ]
      obs_grp1 <- obs_sub[obs_sub$group == grp1, ]
      
      # Arrow target: upper ribbon edge at ~75 min
      idx_ribbon <- which.min(abs(traj_grp1$minute - 75))
      ribbon_x <- traj_grp1$minute[idx_ribbon]
      ribbon_y <- traj_grp1$hi[idx_ribbon]
      
      # Arrow target: top of error bar at 90 min
      obs_90 <- obs_grp1[which.min(abs(obs_grp1$minute - 90)), ]
      errbar_x <- obs_90$minute
      errbar_y <- obs_90$hi
      
      # Positioning reference
      y_max <- max(obs_sub$hi, na.rm = TRUE)
      y_range <- y_max - min(obs_sub$lo, na.rm = TRUE)
      
      # Expand y-axis to make room for labels
      p <- p + ggplot2::coord_cartesian(
        clip = "off",
        ylim = c(NA, y_max + y_range * 0.18)
      )
      
      # Ribbon annotation (left side)
      label_y_ribbon <- y_max + y_range * 0.14
      p <- p +
        ggplot2::annotate(
          "text", x = 5, y = label_y_ribbon,
          label = "True-score 95% CrI",
          size = 3, hjust = 0, fontface = "italic", color = "grey40"
        ) +
        ggplot2::annotate(
          "segment",
          x = 50, y = label_y_ribbon - y_range * 0.02,
          xend = ribbon_x, yend = ribbon_y + y_range * 0.01,
          color = "grey40", linewidth = 0.4,
          arrow = ggplot2::arrow(length = ggplot2::unit(0.12, "cm"), type = "closed")
        )
      
      # Error bar annotation (right side)
      label_y_obs <- y_max + y_range * 0.14
      p <- p +
        ggplot2::annotate(
          "text", x = 120, y = label_y_obs,
          label = "Observed mean \u00b1 95% CI",
          size = 3, hjust = 1, fontface = "italic", color = "grey40"
        ) +
        ggplot2::annotate(
          "segment",
          x = 100, y = label_y_obs - y_range * 0.02,
          xend = errbar_x + 2, yend = errbar_y + y_range * 0.01,
          color = "grey40", linewidth = 0.4,
          arrow = ggplot2::arrow(length = ggplot2::unit(0.12, "cm"), type = "closed")
        )
    }
    
    # Add inset if treat_results available for this analyte
    if (!is.null(treat_results) && spec$auc_label %in% treat_results$index) {
      inset <- plot_treatment_inset(treat_results, spec$auc_label)
      p <- p + patchwork::inset_element(
        inset,
        left = inset_position[1], bottom = inset_position[2],
        right = inset_position[3], top = inset_position[4]
      )
    }
    
    panels[[spec_name]] <- p
  }
  
  # Combine panels with shared legend
  combined <- Reduce(`+`, panels) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  
  combined
}