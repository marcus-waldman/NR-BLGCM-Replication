#' Build per-chain init lists for brms/cmdstanr from a fitted brms model.
#'
#' Samples one random posterior iteration per chain from `fit` and reshapes the
#' draws into the list-of-named-lists form that cmdstanr's `init` argument
#' expects: values on the *constrained* scale, keyed by **Stan parameter
#' block** names with their declared array shapes.
#'
#' Two filtering steps are applied:
#'
#'   1. **Stan parameters block.** We parse `brms::stancode(fit)` to recover
#'      the names declared in the `parameters {}` block (e.g. `b_logglucose0`,
#'      `sd_1`, `z_1`, `L_1`, `Ymi_*`, `Intercept_sigma_*`). brms-renamed
#'      flat names (`b_logglucose0_leg_poly_0_ctrl`,
#'      `r_PTID__logglucose0[14,leg_poly_0]`) are NOT what Stan's init reader
#'      expects, so we use Stan-block names instead.
#'
#'   2. **Saved draws.** brms by default does NOT save the non-centered
#'      standardized random-effects matrix `z_1` or the Cholesky factor `L_1`
#'      (memory savings; the derived `r_1_*` and `cor_*` are saved instead).
#'      Anything declared in the parameters block but not present in
#'      `rstan::extract()` is dropped from the init list — Stan falls back to
#'      its default random init for those, which is fine for non-centered
#'      params naturally near zero.
#'
#' Works whether `fit$fit` is a live `CmdStanMCMC` or has been converted to an
#' rstan `stanfit` (brms's default `read_csv_as_stanfit = TRUE` converts
#' cmdstanr CSV output to a stanfit at fit time so the object survives
#' `saveRDS`/`read_rds`).
#'
#' @param fit A fitted `brmsfit` (any backend).
#' @param n_chains Integer, number of chains the new fit will use (default 4).
#' @param target_vars Optional character vector of Stan-block parameter names
#'   declared by the target model. If supplied, init values are further
#'   filtered to this set. Useful when the target model's parameter block
#'   differs from `fit`'s (e.g. independence variant drops some parameters).
#' @param seed Integer RNG seed for reproducible draw selection (default 42).
#'
#' @return A list of length `n_chains`, each element a named list of parameter
#'   arrays keyed by Stan-block parameter names with declared shapes. Pass
#'   directly as `init = ...` to `brms::brm()`.
make_inits_from_brms <- function(fit, n_chains = 4,
                                 target_vars = NULL, seed = 42) {
  stopifnot(inherits(fit, "brmsfit"))

  set.seed(seed)

  # --- (1) Parse Stan code to find parameters-block declarations ---
  stan_param_vars <- .parse_stan_parameter_names(brms::stancode(fit))
  if (length(stan_param_vars) == 0L) {
    stop("Could not parse any parameters-block declarations from ",
         "brms::stancode(fit).")
  }

  # --- (2) Extract draws keyed by Stan-block names ---
  ext <- rstan::extract(fit$fit, permuted = TRUE)
  keep_names <- intersect(stan_param_vars, names(ext))
  if (length(keep_names) == 0L) {
    stop("No overlap between parameters-block declarations and saved draws. ",
         "Were group-level parameters excluded from saving?")
  }
  ext <- ext[keep_names]

  n_iter   <- dim(ext[[1]])[1]
  draw_idx <- sample.int(n_iter, n_chains, replace = (n_iter < n_chains))

  slice_iter <- function(arr, i) {
    d <- dim(arr)
    if (is.null(d) || length(d) == 1L) return(unname(arr[i]))
    rest_size <- prod(d[-1])
    flat <- seq(from = i, by = d[1], length.out = rest_size)
    array(arr[flat], dim = d[-1])
  }

  lapply(seq_len(n_chains), function(c) {
    i    <- draw_idx[c]
    vals <- lapply(ext, slice_iter, i = i)
    if (!is.null(target_vars)) {
      vals <- vals[intersect(names(vals), target_vars)]
    }
    vals
  })
}


#' Parse the names declared in a Stan `parameters {}` block.
#'
#' Returns the identifier appearing immediately before each declaration's
#' trailing semicolon. Handles type qualifiers like `vector[N]`,
#' `matrix<lower=0>[M, N]`, `array[K] cholesky_factor_corr[M]`, etc.
#'
#' @keywords internal
.parse_stan_parameter_names <- function(stancode) {
  lines <- strsplit(stancode, "\n")[[1]]
  starts <- grep("^\\s*parameters\\s*\\{\\s*$", lines)
  if (length(starts) == 0L) return(character(0))

  start <- starts[1L]
  tail  <- lines[(start + 1L):length(lines)]
  ends  <- which(grepl("^\\}\\s*$", tail))
  if (length(ends) == 0L) return(character(0))
  block <- tail[seq_len(ends[1L] - 1L)]

  block <- gsub("//.*$", "", block)         # strip // line comments
  vars  <- character(0)
  for (ln in block) {
    ln <- trimws(ln)
    if (nchar(ln) == 0L) next
    m <- regmatches(
      ln,
      regexec("([A-Za-z_][A-Za-z0-9_]*)\\s*;\\s*$", ln)
    )[[1L]]
    if (length(m) >= 2L) vars <- c(vars, m[2L])
  }
  unique(vars)
}
