#' Step 5: Figure 6, prior power-scaling robustness of the AUC treatment
#' effects (priorsense; Kallioinen et al., 2024).
#'
#' From the primary fit: (a) the pointwise log-likelihood restricted to
#' observed cells (the mi() rows with a missing response are non-finite in
#' brms::log_lik and are dropped), (b) the per-draw true-score AUC treatment
#' effects (the estimands the figure is about), and (c) the stored lprior.
#' These are assembled into a draws object, cached at
#' derived/priorsense_xdraws.rds, then power-scaled. The figure keeps only the
#' prior-scaling component of the density plot.
#'
#' Run from the repository root:  Rscript scripts/05_make_figure6_priorsense.R
#' Requires the primary fit (step 2). The cache is reused on reruns.

local({
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a) == 1) {
    setwd(normalizePath(file.path(dirname(sub("^--file=", "", a)), "..")))
  }
})
suppressMessages({
  library(brms)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(priorsense)
  library(ggplot2)
})
source("R/_setup.R")

cache = file.path(DERIVED_DIR, "priorsense_xdraws.rds")
quants <- c("Glucose AUC treat effect", "Insulin AUC treat effect")

if (file.exists(cache)) {
  x_draws <- readRDS(cache)
  cat("Loaded cached draws object:", cache, "\n")
} else {
  if (!file.exists(fit_path("primary"))) {
    stop("Missing primary fit at ", fit_path("primary"), ". Run scripts/02_fit_models.R first.")
  }
  fit <- readRDS(fit_path("primary"))

  ll_raw <- brms::log_lik(fit)
  ll_clean <- ll_raw[, is.finite(colSums(ll_raw)), drop = FALSE]
  lp <- posterior::as_draws_df(fit)[["lprior"]]

  nd <- build_clinical_newdata(fit)
  ptid_info <- dplyr::distinct(nd, PTID, treat) |> dplyr::arrange(PTID)
  n_indiv = nrow(ptid_info)
  auc_glu_0 <- true_score_auc_one_resp(fit, nd, "logglucose0", n_indiv)
  auc_ins_0 <- true_score_auc_one_resp(fit, nd, "loginsulin0", n_indiv)
  auc_glu_1 <- true_score_auc_one_resp(fit, nd, "logglucose1", n_indiv)
  auc_ins_1 <- true_score_auc_one_resp(fit, nd, "loginsulin1", n_indiv)
  tau_glucose <- ancova_treat_per_draw(AUC_glucose_1 ~ treat + AUC_glucose_0 + AUC_insulin_0,
                                       auc_glu_0, auc_ins_0, auc_glu_1, auc_ins_1, ptid_info$treat)
  tau_insulin <- ancova_treat_per_draw(AUC_insulin_1 ~ treat + AUC_glucose_0 + AUC_insulin_0,
                                       auc_glu_0, auc_ins_0, auc_glu_1, auc_ins_1, ptid_info$treat)

  ll_df <- as.data.frame(ll_clean)
  names(ll_df) <- paste0("log_lik[", seq_len(ncol(ll_clean)), "]")
  x_draws <- posterior::as_draws_df(data.frame(
    `Glucose AUC treat effect` = tau_glucose,
    `Insulin AUC treat effect` = tau_insulin,
    lprior = lp, ll_df, check.names = FALSE
  ))
  saveRDS(x_draws, cache)
  cat("Computed and cached draws object:", cache, "\n")
}

#-------------------------------------------------------------------------------
# Power-scaling sequence, then keep only the prior-scaling facet of the
# density plot.
#-------------------------------------------------------------------------------
pss <- priorsense::powerscale_sequence(x_draws)
p_full <- priorsense::powerscale_plot_dens(pss, variable = quants, help_text = FALSE)

keep_prior <- function(df) {
  return(df[as.character(df$component) == "prior", , drop = FALSE])
}
p_prior <- p_full
p_prior$data <- keep_prior(p_prior$data)
for (i in seq_along(p_prior$layers)) {
  ld <- p_prior$layers[[i]]$data
  if (is.data.frame(ld) && "component" %in% names(ld)) {
    p_prior$layers[[i]]$data <- keep_prior(ld)
  }
}
p_prior <- p_prior +
  ggplot2::facet_wrap(~ variable, scales = "free", nrow = 1) +
  ggplot2::labs(
    title = "Robustness of treatment-effect estimates to prior power-scaling",
    subtitle = "Posterior densities at prior scaled 0.8x / 1x / 1.25x; overlap = prior-insensitive",
    x = NULL, y = "Posterior density"
  ) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(file.path(FIG_DIR, "Fig6_priorsense_prior_robustness.png"),
                p_prior, width = 9, height = 4.5, dpi = 300)
saveRDS(p_prior, file.path(FIG_DIR, "Fig6_priorsense_prior_robustness.rds"))

# Console companion: priorsense's sensitivity diagnostics for the two estimands.
cat("\n== priorsense::powerscale_sensitivity ==\n")
print(priorsense::powerscale_sensitivity(x_draws, variable = quants))
cat("\nWrote", file.path(FIG_DIR, "Fig6_priorsense_prior_robustness.png"), "\n")
