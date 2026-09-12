#' Fit Legendre Polynomials to External OGTT Reference Data (Carlson 2025, v2)
#'
#' Audit / diagnostic script. The runtime priors used by the model are derived
#' inside `utils/prior_specification.R::derive_legendre_priors_from_csv()` from
#' the SAME v2 CSV consumed here, using the SAME method. This script exists to
#' make the derivation inspectable: it writes a coefficient table and a fit
#' plot to disk, and prints the priors that would be emitted at runtime.
#'
#' Source data: `external_ogtt_glucose_reference v2.csv`
#'   - Wide format with two independent visual-rater readings per group:
#'     `gdm_my_estimate`, `gdm_other_estimate`,
#'     `no_gdm_my_estimate`, `no_gdm_other_estimate` (plus precomputed `*_avg`).
#'   - GDM N=15, No-GDM N=29 from the published study figure.
#'
#' Derivation:
#'   1. Restrict to t in [0, 120] min; map to Legendre domain via t/60 - 1.
#'   2. log-transform glucose values.
#'   3. For each group (GDM, No-GDM), stack both raters' rows and fit
#'      `lm(log_value ~ P0 + P1 + P2 + P3 + P4 - 1)`. The lm SE per coefficient
#'      therefore reflects digitization / rater noise.
#'   4. prior_mean = coef_gdm (GDM-centered).
#'      prior_sd   = sqrt(se_gdm^2 + ((coef_gdm - coef_no_gdm)/1.5)^2)
#'      (rater noise + between-group spread combined in quadrature).
#'
#' P0 is reported here for completeness but is NOT used by the runtime
#' Legendre prior (P0 receives a data-scaled prior in
#' `create_joined_growth_curve_priors_n()`).

library(dplyr)
library(tidyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# 0. Resolve script directory so this runs whether sourced or Rscript'd.
# -----------------------------------------------------------------------------
script_dir <- tryCatch({
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    # Fall back to commandArgs heuristic for Rscript
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) == 1L) {
      normalizePath(dirname(sub("^--file=", "", file_arg)))
    } else {
      getwd()
    }
  }
}, error = function(e) getwd())
setwd(script_dir)

csv_path <- "external_ogtt_glucose_reference v2.csv"

cat("=======================================================================\n")
cat("FITTING LEGENDRE POLYNOMIALS TO CARLSON 2025 OGTT GLUCOSE DATA (v2)\n")
cat("=======================================================================\n\n")
cat("Source CSV:", csv_path, "\n\n")

# -----------------------------------------------------------------------------
# 1. Load v2 CSV and reshape to long per-rater format.
# -----------------------------------------------------------------------------
raw <- read.csv(csv_path, check.names = TRUE, stringsAsFactors = FALSE)

needed <- c("time_min",
            "gdm_my_estimate", "gdm_other_estimate",
            "no_gdm_my_estimate", "no_gdm_other_estimate")
missing_cols <- setdiff(needed, names(raw))
if (length(missing_cols) > 0) {
  stop("v2 CSV missing columns: ", paste(missing_cols, collapse = ", "))
}

long <- bind_rows(
  tibble(time_min = raw$time_min, rater = "rater_1",
         group = "gdm",    value = raw$gdm_my_estimate),
  tibble(time_min = raw$time_min, rater = "rater_2",
         group = "gdm",    value = raw$gdm_other_estimate),
  tibble(time_min = raw$time_min, rater = "rater_1",
         group = "no_gdm", value = raw$no_gdm_my_estimate),
  tibble(time_min = raw$time_min, rater = "rater_2",
         group = "no_gdm", value = raw$no_gdm_other_estimate)
) %>%
  filter(time_min >= 0, time_min <= 120) %>%
  mutate(
    t_scaled  = time_min / 60 - 1,
    log_value = log(value)
  )

cat("Stacked-rater long table (head):\n")
print(head(long))
cat("\nTimepoints used:",
    paste(sort(unique(long$time_min)), collapse = ", "), "min\n")
cat("Raters:", paste(sort(unique(long$rater)), collapse = ", "), "\n")
cat("Groups:", paste(sort(unique(long$group)), collapse = ", "), "\n")
cat("Rows per group:",
    paste(sapply(split(long, long$group), nrow), collapse = " / "), "\n\n")

# -----------------------------------------------------------------------------
# 2. Legendre basis P0..P4.
# -----------------------------------------------------------------------------
legendre_basis <- function(x) {
  data.frame(
    P0 = rep(1, length(x)),
    P1 = x,
    P2 = (3 * x^2 - 1) / 2,
    P3 = (5 * x^3 - 3 * x) / 2,
    P4 = (35 * x^4 - 30 * x^2 + 3) / 8
  )
}
long <- cbind(long, legendre_basis(long$t_scaled))

# -----------------------------------------------------------------------------
# 3. Fit per group with stacked rater rows.
# -----------------------------------------------------------------------------
fit_group <- function(g) {
  dat <- long %>% filter(group == g)
  fit <- lm(log_value ~ P0 + P1 + P2 + P3 + P4 - 1, data = dat)
  smry <- summary(fit)$coefficients
  tibble(
    term = rownames(smry),
    coef = smry[, "Estimate"],
    se   = smry[, "Std. Error"]
  )
}

gdm_fit    <- fit_group("gdm")
no_gdm_fit <- fit_group("no_gdm")

cat("--- GDM fit (both raters stacked) ---\n")
print(gdm_fit, n = 5)
cat("\n--- No-GDM fit (both raters stacked) ---\n")
print(no_gdm_fit, n = 5)
cat("\n")

# -----------------------------------------------------------------------------
# 4. Build prior table: GDM-centered, quadrature SD.
# -----------------------------------------------------------------------------
prior_tbl <- gdm_fit %>%
  rename(coef_gdm = coef, se_gdm = se) %>%
  inner_join(
    no_gdm_fit %>% rename(coef_no_gdm = coef, se_no_gdm = se),
    by = "term"
  ) %>%
  mutate(
    group_spread = abs(coef_gdm - coef_no_gdm) / 1.5,
    prior_mean   = coef_gdm,
    prior_sd     = sqrt(se_gdm^2 + group_spread^2)
  )

cat("=======================================================================\n")
cat("LEGENDRE COEFFICIENT TABLE (log-glucose scale)\n")
cat("  prior_mean = coef_gdm\n")
cat("  prior_sd   = sqrt(se_gdm^2 + ((coef_gdm - coef_no_gdm)/1.5)^2)\n")
cat("=======================================================================\n\n")
print(prior_tbl, n = 5, width = 120)
cat("\n")

write.csv(prior_tbl,
          file = "legendre_coefficients_external_ogtt.csv",
          row.names = FALSE)
cat("Wrote: legendre_coefficients_external_ogtt.csv\n")

# Supplement-ready S2 Table 2: P1..P4 only, 2 decimal places, columns matching
# the manuscript table (Term, GDM, No GDM, Prior Mean, Prior SD).
s2_table2 <- prior_tbl %>%
  filter(term != "P0") %>%
  transmute(
    Term       = term,
    GDM        = sprintf("%.2f", coef_gdm),
    `No GDM`   = sprintf("%.2f", coef_no_gdm),
    `Prior Mean` = sprintf("%.2f", prior_mean),
    `Prior SD` = sprintf("%.2f", prior_sd)
  )
write.csv(s2_table2,
          file = "S2_table2_glucose_legendre_priors.csv",
          row.names = FALSE)
cat("Wrote: S2_table2_glucose_legendre_priors.csv\n\n")

# -----------------------------------------------------------------------------
# 5. Print the priors that runtime would emit (P1..P4 only; P0 is data-scaled).
# -----------------------------------------------------------------------------
cat("=======================================================================\n")
cat("PRIORS EMITTED AT RUNTIME (P1-P4; P0 is data-scaled elsewhere)\n")
cat("=======================================================================\n\n")
for (i in seq_len(nrow(prior_tbl))) {
  tm <- prior_tbl$term[i]
  if (tm == "P0") next
  cat(sprintf("%s: Normal(%.4f, %.4f)\n",
              tm,
              prior_tbl$prior_mean[i],
              prior_tbl$prior_sd[i]))
}
cat("\n")

# -----------------------------------------------------------------------------
# 6. Fit-plot with prior-SD spaghetti.
#    Draws from N(prior_mean, prior_sd) for P1..P4 (GDM-centered); P0 held at
#    the GDM fitted value. Each draw evaluated on a dense time grid and plotted
#    as a faint line to visualize prior uncertainty in the trajectory shape.
# -----------------------------------------------------------------------------
set.seed(42)
n_draws <- 400
kappa_grid <- c(1, 10, 100)  # Prior-SD scaling factors for the three panels.
t_grid_scaled <- seq(-1, 1, length.out = 200)
basis_grid <- legendre_basis(t_grid_scaled)
basis_mat  <- as.matrix(basis_grid[, c("P0", "P1", "P2", "P3", "P4")])

# Sample P1..P4 from prior at each kappa; hold P0 at fitted GDM intercept.
prior_for <- function(tm) prior_tbl[prior_tbl$term == tm, , drop = FALSE]
P0_fixed <- prior_for("P0")$coef_gdm
sample_draws <- function(kappa) {
  draws <- matrix(NA_real_, nrow = n_draws, ncol = 5,
                  dimnames = list(NULL, c("P0", "P1", "P2", "P3", "P4")))
  draws[, "P0"] <- P0_fixed
  for (tm in c("P1", "P2", "P3", "P4")) {
    r <- prior_for(tm)
    draws[, tm] <- rnorm(n_draws, mean = r$prior_mean, sd = r$prior_sd * kappa)
  }
  draws
}

# Group geometric means (raw scale) used to center each curve.
P0_gdm    <- coef(lm(log_value ~ P0 + P1 + P2 + P3 + P4 - 1,
                     data = long %>% filter(group == "gdm")))["P0"]
P0_no_gdm <- coef(lm(log_value ~ P0 + P1 + P2 + P3 + P4 - 1,
                     data = long %>% filter(group == "no_gdm")))["P0"]
mean_gdm    <- exp(P0_gdm)
mean_no_gdm <- exp(P0_no_gdm)

# Spaghetti per kappa: GDM-centered prior on P1..P4, P0 held at GDM intercept,
# then centered on the GDM mean so y-axis is deviation in mg/dL.
build_spaghetti <- function(kappa) {
  draws <- sample_draws(kappa)
  # Log-scale shape: subtract P0 to remove the intercept; the result is the
  # polynomial deviation in log-glucose units.
  mat <- basis_mat %*% t(draws) - P0_fixed
  sp  <- as.data.frame(t(mat))
  sp$draw_id <- seq_len(n_draws)
  tidyr::pivot_longer(sp, -draw_id, names_to = "grid_idx", values_to = "delta") %>%
    mutate(grid_idx = as.integer(sub("V", "", grid_idx)),
           time_min = (t_grid_scaled[grid_idx] + 1) * 60,
           kappa    = kappa)
}
spaghetti_long <- bind_rows(lapply(kappa_grid, build_spaghetti)) %>%
  mutate(kappa_label = factor(sprintf("κ = %d", kappa),
                              levels = sprintf("κ = %d", kappa_grid)))

# Group-specific fitted curves, each centered on its own group mean.
coef_mat <- list(
  gdm    = setNames(gdm_fit$coef,    gdm_fit$term),
  no_gdm = setNames(no_gdm_fit$coef, no_gdm_fit$term)
)
pred_grid <- bind_rows(
  tibble(time_min = (t_grid_scaled + 1) * 60,
         group = "gdm",
         delta_pred = (basis_mat %*% coef_mat$gdm)[, 1] - P0_gdm),
  tibble(time_min = (t_grid_scaled + 1) * 60,
         group = "no_gdm",
         delta_pred = (basis_mat %*% coef_mat$no_gdm)[, 1] - P0_no_gdm)
)
# Replicate fitted curves across each kappa panel.
pred_grid <- bind_rows(lapply(kappa_grid, function(k) pred_grid %>% mutate(kappa = k))) %>%
  mutate(kappa_label = factor(sprintf("κ = %d", kappa),
                              levels = sprintf("κ = %d", kappa_grid)))

# Observed points, centered on their group mean and replicated across panels.
long_centered <- long %>%
  mutate(delta = log_value - if_else(group == "gdm", as.numeric(P0_gdm), as.numeric(P0_no_gdm)))
long_centered <- bind_rows(lapply(kappa_grid, function(k) long_centered %>% mutate(kappa = k))) %>%
  mutate(kappa_label = factor(sprintf("κ = %d", kappa),
                              levels = sprintf("κ = %d", kappa_grid)))

# Okabe-Ito colorblind-safe palette: vermillion + blue.
ok_gdm    <- "#D55E00"
ok_no_gdm <- "#0072B2"

p <- ggplot() +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
  # Prior-SD spaghetti (GDM-centered, mean-subtracted).
  geom_line(data = spaghetti_long,
            aes(x = time_min, y = delta, group = draw_id),
            color = ok_gdm, alpha = 0.04, linewidth = 0.25) +
  # Fitted group curves, mean-centered.
  geom_line(data = pred_grid,
            aes(x = time_min, y = delta_pred, color = group),
            linewidth = 1.1) +
  # Per-rater observed points, mean-centered.
  geom_point(data = long_centered,
             aes(x = time_min, y = delta,
                 color = group, shape = rater),
             size = 2.6, alpha = 0.9) +
  scale_color_manual(values = c(gdm = ok_gdm, no_gdm = ok_no_gdm),
                     labels = c(gdm = "GDM", no_gdm = "No-GDM")) +
  scale_shape_manual(values = c(rater_1 = 16, rater_2 = 17)) +
  facet_wrap(~ kappa_label, nrow = 1, scales = "free_y") +
  labs(
    title = "Legendre P0-P4 Fits to Carlson 2025 OGTT Glucose",
    x = "Time (min)",
    y = "log(glucose / group geometric mean)",
    color = "Group", shape = "Rater"
  ) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold"))

ggsave("legendre_fit_external_ogtt.png", p,
       width = 14, height = 5, dpi = 150)
cat("Wrote: legendre_fit_external_ogtt.png\n")

cat("\nScript complete.\n")
