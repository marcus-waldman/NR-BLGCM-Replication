# ============================================================================
# TREATMENT EFFECT INSET
# ============================================================================

#' Plot Compact Treatment Effect Inset
#'
#' Creates a small horizontal pointrange plot for embedding as an inset in
#' trajectory panels. Shows observed vs true-score AUC treatment effects
#' with numeric annotations.
#'
#' @param treat_results data.frame with columns: index, method, estimate,
#'   conf.low, conf.high
#' @param auc_label Character, the index name to filter from treat_results
#'   (e.g., "Glucose AUC")
#' @param base_size Numeric base font size (default 8)
#'
#' @return A ggplot object suitable for use with patchwork::inset_element()
#'
#' @export
plot_treatment_inset <- function(treat_results, auc_label,
                                 base_size = 8) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required.")
  }

  inset_data <- treat_results[treat_results$index == auc_label, ]

  # Simplify method labels for inset display
  inset_data$method_short <- ifelse(
    grepl("True-score", inset_data$method), "True score", "Observed"
  )
  inset_data$method_short <- factor(inset_data$method_short,
                                    levels = c("True score", "Observed"))

  # Format annotation text: "est [lo, hi]"
  inset_data$label <- sprintf("%.1f [%.1f, %.1f]",
                              inset_data$estimate,
                              inset_data$conf.low,
                              inset_data$conf.high)

  ggplot2::ggplot(inset_data, ggplot2::aes(y = method_short, x = estimate)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40",
                        linewidth = 0.4) +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = conf.low, xmax = conf.high),
      size = 0.4, linewidth = 0.7, colour = "black"
    ) +
    ggplot2::geom_label(
      ggplot2::aes(label = label),
      vjust = -0.8, size = 2.8, show.legend = FALSE,
      fill = "white", colour = "black",
      linewidth = 0, label.padding = ggplot2::unit(1, "pt")
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = "AUC Treatment Effect") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = "black",
                                              linewidth = 0.7),
      plot.title = ggplot2::element_text(size = 8, hjust = 0.5, face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7.5, lineheight = 0.9, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 6),
      plot.margin = ggplot2::margin(3, 5, 3, 3)
    )
}
