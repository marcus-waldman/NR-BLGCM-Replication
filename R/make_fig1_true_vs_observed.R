#' Figure 1: Conceptual figure -- True Score vs Observed Trajectory
#'
#' One observed measurement per timepoint with SEM error bars on the
#' log-glucose scale; dual y-axis (log glucose primary, raw mg/dL secondary).
#' Posterior measurement-error estimate (~4% CV) matches the Bayesian
#' compromise illustrated in Figure 2.
#'
#' Pure conceptual figure: no dependence on fitted models.
#'
#' Color palette: JAMA (ggsci) -- blue = knowledge / true score / prior,
#' orange = data / observations / likelihood, dark teal = posterior / SEM.
#' Shared with Figure 2.
make_fig1_true_vs_observed <- function(
    output_png = "revisions/figures/Fig1_true_vs_observed.png",
    output_rds = "revisions/figures/Fig1_true_vs_observed.rds",
    seed = 42
) {

  set.seed(seed)

  col_knowledge = "#00A1D5"
  col_data      = "#DF8F44"
  col_synthesis = "#374E55"

  time_fine <- seq(0, 120, by = 1)

  true_trajectory_raw <- function(t) {
    return(90 + 80 * exp(-0.02 * (t - 45)^2 / 100) * (1 - exp(-t / 10)))
  }
  true_trajectory_log <- function(t) {
    return(log(true_trajectory_raw(t)))
  }

  time_clinical = c(0, 30, 60, 90, 120)
  true_log <- true_trajectory_log(time_clinical)
  true_raw <- true_trajectory_raw(time_clinical)

  # Posterior measurement-error SD on log scale (~4% CV); see Fig 2.
  sigma_log = 0.040

  # Hand-picked offsets for visual clarity (all within ~1.5 posterior SEM).
  obs_log <- true_log + c(0.03, -0.05, 0.04, -0.025, 0.015)
  obs_raw <- exp(obs_log)

  sem_upper_log <- obs_log + sigma_log
  sem_lower_log <- obs_log - sigma_log

  smooth_data <- data.frame(
    time = time_fine,
    log_glucose = true_trajectory_log(time_fine)
  )

  point_data <- data.frame(
    time = time_clinical,
    true_log = true_log,
    obs_log = obs_log,
    sem_upper = sem_upper_log,
    sem_lower = sem_lower_log
  )

  raw_breaks = c(80, 90, 100, 120, 140, 160, 180)
  log_breaks <- log(raw_breaks)

  p <- ggplot2::ggplot() +
    ggplot2::geom_errorbar(
      data = point_data,
      ggplot2::aes(x = time, ymin = sem_lower, ymax = sem_upper),
      color = col_synthesis, width = 4, linewidth = 0.6, alpha = 0.5
    ) +
    ggplot2::geom_line(
      data = smooth_data,
      ggplot2::aes(x = time, y = log_glucose),
      color = col_knowledge, linewidth = 1.2
    ) +
    ggplot2::geom_segment(
      data = point_data,
      ggplot2::aes(x = time, xend = time, y = true_log, yend = obs_log),
      color = "grey50", linewidth = 0.4, linetype = "dotted"
    ) +
    ggplot2::geom_point(
      data = point_data,
      ggplot2::aes(x = time, y = true_log),
      color = col_knowledge, size = 3, shape = 15
    ) +
    ggplot2::geom_point(
      data = point_data,
      ggplot2::aes(x = time, y = obs_log),
      color = col_data, size = 3, shape = 16
    ) +
    ggplot2::annotate(
      "text", x = 108, y = point_data$true_log[4] + 0.06,
      label = expression("true score," ~ eta),
      color = col_knowledge, hjust = 0, size = 3.5, fontface = "italic"
    ) +
    ggplot2::annotate(
      "segment",
      x = 107, y = point_data$true_log[4] + 0.055,
      xend = 91.5, yend = point_data$true_log[4] + 0.005,
      color = col_knowledge, linewidth = 0.4,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.12, "cm"), type = "closed")
    ) +
    ggplot2::annotate(
      "text", x = 108, y = point_data$obs_log[4] - 0.06,
      label = "observed measurement",
      color = col_data, hjust = 0, size = 3.5, fontface = "italic"
    ) +
    ggplot2::annotate(
      "segment",
      x = 107, y = point_data$obs_log[4] - 0.055,
      xend = 91.5, yend = point_data$obs_log[4] - 0.005,
      color = col_data, linewidth = 0.4,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.12, "cm"), type = "closed")
    ) +
    ggplot2::annotate(
      "text", x = 68, y = sem_lower_log[4] - 0.04,
      label = expression("± SEM (~4% CV)"),
      color = col_synthesis, hjust = 0.5, size = 3.5, fontface = "italic",
      alpha = 0.5
    ) +
    ggplot2::annotate(
      "segment",
      x = 78, y = sem_lower_log[4] - 0.035,
      xend = 89, yend = sem_lower_log[4] + 0.005,
      color = col_synthesis, linewidth = 0.4, alpha = 0.5,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.12, "cm"), type = "closed")
    ) +
    ggplot2::scale_y_continuous(
      name = "Log glucose",
      sec.axis = ggplot2::sec_axis(
        transform = ~ exp(.),
        name = "Glucose (mg/dL)",
        breaks = raw_breaks
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = time_clinical,
      labels = ifelse(time_clinical == 0, "Fasting\n(0 min)", paste0(time_clinical, " min")),
      expand = ggplot2::expansion(mult = c(0.02, 0.25))
    ) +
    ggplot2::labs(x = "Time after glucose load") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.y.left = ggplot2::element_text(color = "grey30"),
      axis.title.y.right = ggplot2::element_text(color = "grey30"),
      plot.margin = ggplot2::margin(10, 20, 10, 10)
    )

  if (!is.null(output_png)) {
    ggplot2::ggsave(output_png, plot = p, width = 8, height = 4.5, dpi = 300)
  }
  if (!is.null(output_rds)) {
    saveRDS(p, output_rds)
  }

  return(p)
}
