# ============================================================================
# TRAJECTORY PANEL (SINGLE ANALYTE)
# ============================================================================

#' Plot Single-Analyte Trajectory Panel
#'
#' Creates a single panel with true-score ribbon + mean line overlaid with
#' observed pointrange at clinical timepoints.
#'
#' @param traj_data data.frame of true-score trajectories (from compute_true_score_trajectories)
#'   filtered to one analyte. Columns: minute, mean, lo, hi, group, analyte.
#' @param obs_data data.frame of observed group means (from compute_observed_group_means)
#'   filtered to one analyte. Same columns.
#' @param y_label Character y-axis label (e.g., "Glucose (mg/dL)")
#' @param group_colours Named character vector of group colours
#' @param group_shapes Named integer vector of group point shapes
#' @param base_size Numeric base font size (default 12)
#'
#' @return A ggplot object
#'
#' @export
plot_trajectory_panel <- function(traj_data, obs_data, y_label,
                                  group_colours = c("LC/Conventional" = "#DF8F44",
                                                    "CHOICE" = "#00A1D5"),
                                  group_shapes = c("LC/Conventional" = 16,
                                                   "CHOICE" = 17),
                                  base_size = 12) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required.")
  }
  
  ggplot2::ggplot() +
    # True-score ribbon
    ggplot2::geom_ribbon(
      data = traj_data,
      ggplot2::aes(x = minute, ymin = lo, ymax = hi, fill = group),
      alpha = 0.15
    ) +
    # Ribbon edge lines (95% CrI bounds)
    ggplot2::geom_line(
      data = traj_data,
      ggplot2::aes(x = minute, y = lo, colour = group),
      linewidth = 0.3, alpha = 0.5
    ) +
    ggplot2::geom_line(
      data = traj_data,
      ggplot2::aes(x = minute, y = hi, colour = group),
      linewidth = 0.3, alpha = 0.5
    ) +
    # True-score smooth line (posterior mean)
    ggplot2::geom_line(
      data = traj_data,
      ggplot2::aes(x = minute, y = mean, colour = group),
      linewidth = 0.9
    ) +
    # Observed means with error bars + distinct shapes
    ggplot2::geom_pointrange(
      data = obs_data,
      ggplot2::aes(x = minute, y = mean, ymin = lo, ymax = hi,
                   colour = group, shape = group),
      position = ggplot2::position_dodge(width = 4),
      size = 0.4, linewidth = 0.6
    ) +
    ggplot2::scale_colour_manual(values = group_colours) +
    ggplot2::scale_fill_manual(values = group_colours) +
    ggplot2::scale_shape_manual(values = group_shapes) +
    ggplot2::labs(
      x = "Minutes",
      y = y_label,
      colour = "Group",
      fill = "Group",
      shape = "Group"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom")
}
