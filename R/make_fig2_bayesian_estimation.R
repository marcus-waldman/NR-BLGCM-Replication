#' Figure 2: Conceptual figure -- Prior / Likelihood / Posterior
#'
#' Shows how Bayesian estimation combines prior knowledge with data, using
#' the measurement-error SD (log glucose scale) as the running example:
#'   Prior:      N(0.030, 0.012) -- biological-variation literature (~3% CV)
#'   Likelihood: N(0.062, 0.018) -- study data (~6% CV)
#'   Posterior:  precision-weighted compromise (~4% CV)
#'
#' Pure conceptual figure: no dependence on fitted models.
#'
#' Color palette: JAMA -- blue = knowledge / prior, orange = data / likelihood,
#' dark teal = posterior. Shared with Figure 1.
make_fig2_bayesian_estimation <- function(
    output_png = "revisions/figures/Fig2_bayesian_estimation.png",
    output_rds = "revisions/figures/Fig2_bayesian_estimation.rds",
    seed = 42
) {

  set.seed(seed)

  col_knowledge = "#00A1D5"
  col_data      = "#DF8F44"
  col_synthesis = "#374E55"

  x <- seq(0.005, 0.12, length.out = 500)

  prior_mean = 0.030
  prior_sd   = 0.012
  prior_density <- dnorm(x, prior_mean, prior_sd)

  lik_mean = 0.062
  lik_sd   = 0.018
  lik_density <- dnorm(x, lik_mean, lik_sd)

  post_prec = 1/prior_sd^2 + 1/lik_sd^2
  post_mean = (prior_mean/prior_sd^2 + lik_mean/lik_sd^2) / post_prec
  post_sd   = 1/sqrt(post_prec)
  post_density <- dnorm(x, post_mean, post_sd)

  df <- rbind(
    data.frame(x = x, density = prior_density,
               Distribution = "Prior\n(biological variation literature)"),
    data.frame(x = x, density = lik_density,
               Distribution = "Likelihood\n(study data)"),
    data.frame(x = x, density = post_density,
               Distribution = "Posterior\n(combined estimate)")
  )

  df$Distribution <- factor(
    df$Distribution,
    levels = c(
      "Prior\n(biological variation literature)",
      "Likelihood\n(study data)",
      "Posterior\n(combined estimate)"
    )
  )

  # Posterior legend swatch: midpoint of #374E55 at 50% alpha on white,
  # to match its faded appearance in the plot.
  col_post_legend = "#9BA5A8"
  colors <- c(
    "Prior\n(biological variation literature)" = col_knowledge,
    "Likelihood\n(study data)"                 = col_data,
    "Posterior\n(combined estimate)"           = col_post_legend
  )

  df_full <- df[df$Distribution != "Posterior\n(combined estimate)", ]
  df_post <- df[df$Distribution == "Posterior\n(combined estimate)", ]

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = df_full,
      ggplot2::aes(x = x, ymin = 0, ymax = density, fill = Distribution),
      alpha = 0.12, color = NA
    ) +
    ggplot2::geom_line(
      data = df_full,
      ggplot2::aes(x = x, y = density, color = Distribution),
      linewidth = 1.1
    ) +
    ggplot2::geom_ribbon(
      data = df_post,
      ggplot2::aes(x = x, ymin = 0, ymax = density, fill = Distribution),
      alpha = 0.06, color = NA
    ) +
    ggplot2::geom_line(
      data = df_post,
      ggplot2::aes(x = x, y = density, color = Distribution),
      linewidth = 1.1, alpha = 0.5
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_x_continuous(
      labels = function(x) { return(paste0(round(x * 100, 1), "%")) },
      breaks = seq(0, 0.12, by = 0.02)
    ) +
    ggplot2::labs(
      x = expression("Residual SD (log scale) " %->% " approximate CV"),
      y = "Density",
      color = NULL,
      fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(10, 15, 10, 10)
    ) +
    ggplot2::geom_vline(xintercept = prior_mean, color = col_knowledge, linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = lik_mean,   color = col_data,      linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = post_mean,  color = col_synthesis, linetype = "dashed", linewidth = 0.5, alpha = 0.5)

  if (!is.null(output_png)) {
    ggplot2::ggsave(output_png, plot = p, width = 6, height = 4, dpi = 300)
  }
  if (!is.null(output_rds)) {
    saveRDS(p, output_rds)
  }

  cat(sprintf("Fig 2 posterior mean: %.4f (~%.1f%% CV)\n", post_mean, post_mean * 100))

  return(p)
}
