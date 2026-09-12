#' Forest plot of the cross-adjusted ANCOVA treatment effect across prior
#' SD scales. Single function used by both main.R (saving as Fig5...) and the
#' standalone driver in revisions/figures/, so the plot definition lives in
#' exactly one place.
#'
#' @param prior_sens_df data.frame with columns:
#'   analyte         (chr or factor; "Glucose"/"Insulin" or "Glucose AUC"/"Insulin AUC")
#'   prior_sd_scale  (numeric; reference scale is detected as == 1)
#'   est, ci_lb, ci_ub (numeric; AUC treat coefficient and 95% bounds)
#'   One row per (analyte x prior_sd_scale). Caller is responsible for any
#'   filtering (e.g. dropping the independence model row, prior_sd_scale = 999).
#' @param output_path optional file path. If non-NULL, ggsave saves the plot.
#' @param width,height,dpi ggsave parameters, only used when output_path is set.
#' @return ggplot object (invisibly when output_path is given).
make_prior_sensitivity_forest <- function(prior_sens_df,
                                          output_path = NULL,
                                          width  = 7.5,
                                          height = 4.5,
                                          dpi    = 300) {
  for (pkg in c("ggplot2", "ggtext", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required.")
    }
  }

  # Build label ordering: diffuse priors on top (descending numeric scale).
  scales_desc  <- sort(unique(prior_sens_df$prior_sd_scale), decreasing = TRUE)
  level_labels <- paste0(scales_desc, "×")  # use multiplication sign

  df <- prior_sens_df
  # Normalize analyte to "Glucose AUC"/"Insulin AUC" factor for facet titles.
  if (is.factor(df$analyte) &&
      all(levels(df$analyte) %in% c("Glucose AUC", "Insulin AUC"))) {
    # already in canonical form
  } else {
    df$analyte <- factor(
      ifelse(grepl("AUC", df$analyte), df$analyte, paste(df$analyte, "AUC")),
      levels = c("Glucose AUC", "Insulin AUC")
    )
  }
  df$is_ref      <- df$prior_sd_scale == 1
  df$prior_label <- factor(paste0(df$prior_sd_scale, "×"),
                           levels = level_labels)

  # Okabe-Ito colorblind-safe: vermillion highlights the reference prior;
  # all non-reference rows are softened to grey50.
  col_ref   <- "#D55E00"
  col_other <- "grey50"

  # Per-row "EAP [low, high]" label, rounded to the nearest tenth.
  df$label_txt <- sprintf("%.1f [%.1f, %.1f]", df$est, df$ci_lb, df$ci_ub)

  # Reference rows (analytic model). Numeric y_pos lets us draw the leader
  # arrow + label on the (otherwise discrete) y-axis cleanly per facet.
  ref_df <- df[df$is_ref, , drop = FALSE]
  ref_df$y_pos <- as.numeric(ref_df$prior_label)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = est, y = prior_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lb, xmax = ci_ub, color = is_ref, linewidth = is_ref),
      width = 0.25, orientation = "y"
    ) +
    ggplot2::geom_point(ggplot2::aes(color = is_ref, size = is_ref)) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = col_other, `TRUE` = col_ref), guide = "none"
    ) +
    ggplot2::scale_linewidth_manual(
      values = c(`FALSE` = 0.6, `TRUE` = 1.2), guide = "none"
    ) +
    ggplot2::scale_size_manual(
      values = c(`FALSE` = 2.6, `TRUE` = 4.0), guide = "none"
    ) +
    ggplot2::facet_wrap(~ analyte, scales = "free_x") +
    ggplot2::labs(
      title = "Prior Sensitivity Analysis",
      x     = "Adjusted Treatment Effect",
      y     = "highly diffuse  ←——  **Prior SD scaling**  ——→  highly informative"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      strip.text       = ggplot2::element_text(face = "bold"),
      axis.title.y     = ggtext::element_markdown(size = 10, angle = 90),
      axis.text.y      = ggplot2::element_text(
        color = ifelse(levels(df$prior_label) == "1×", col_ref, col_other),
        face  = ifelse(levels(df$prior_label) == "1×", "bold",   "plain")
      ),
      plot.margin      = ggplot2::margin(t = 8, r = 10, b = 8, l = 8)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  # Per-bar "EAP [low, high]" labels (all rows). Vermillion for the analytic
  # model, grey50 for the sensitivity rungs. Placed just above each point.
  p <- p +
    ggplot2::geom_text(
      data = df,
      ggplot2::aes(
        x     = est,
        y     = as.numeric(prior_label) + 0.22,
        label = label_txt,
        color = is_ref
      ),
      size        = 3.4,
      inherit.aes = FALSE,
      show.legend = FALSE
    )

  if (nrow(ref_df) > 0) {
    # Bumped higher than before so the leader+label clears the EAP text above
    # the 1x bar (which sits at y_pos + 0.22).
    p <- p +
      ggplot2::geom_segment(
        data = ref_df,
        ggplot2::aes(x = est, xend = est,
                     y = y_pos + 0.68, yend = y_pos + 0.45),
        inherit.aes = FALSE,
        color       = col_ref,
        linewidth   = 0.5,
        arrow       = grid::arrow(length = grid::unit(0.07, "inches"),
                                  ends = "last", type = "closed")
      ) +
      ggplot2::geom_text(
        data = ref_df,
        ggplot2::aes(x = est, y = y_pos + 0.86, label = "Analytic Model"),
        inherit.aes = FALSE,
        color       = col_ref,
        fontface    = "bold",
        size        = 3.4
      )
  }

  if (!is.null(output_path)) {
    ggplot2::ggsave(filename = output_path, plot = p,
                    width = width, height = height, dpi = dpi)
    return(invisible(p))
  }
  p
}
