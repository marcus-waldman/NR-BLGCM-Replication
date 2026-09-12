#' Evaluate Legendre Polynomial of Degree k
#'
#' Computes the k-th degree Legendre polynomial P_k(t) for t ∈ [-1, 1].
#' Uses the orthopolynom package to generate polynomial coefficients via
#' Bonnet's recursion formula.
#'
#' @param t Numeric vector of values in [-1, 1]
#' @param k Integer degree of polynomial (0, 1, 2, 3, 4, ...)
#'
#' @return Numeric vector of P_k(t) values
#'
#' @examples
#' t <- seq(-1, 1, length.out = 100)
#' p1 <- leg_poly(t, 1)  # Linear: P_1(t) = t
#' p2 <- leg_poly(t, 2)  # Quadratic: P_2(t) = (3t² - 1)/2
#'
#' @export
leg_poly <- function(t, k) {
  if (!requireNamespace("orthopolynom", quietly = TRUE)) {
    stop("Package 'orthopolynom' required but not installed.\n",
         "Install with: install.packages('orthopolynom')")
  }
  
  if (k < 0 || k != round(k)) {
    stop("k must be a non-negative integer")
  }
  
  if (any(t < -1 | t > 1, na.rm = TRUE)) {
    warning("Some t values are outside [-1, 1]. Legendre polynomials are defined on [-1, 1].")
  }
  
  # Generate Legendre polynomial coefficients up to degree k
  # Returns list of polynomial.functions objects
  leg_list <- orthopolynom::legendre.polynomials(n = k, normalized = FALSE)
  
  # Extract the k-th polynomial (list is 0-indexed: [[1]] = P_0, [[2]] = P_1, etc.)
  poly_k <- leg_list[[k + 1]]
  print(c("poly_k",poly_k)) #sanity check
  
  # Evaluate polynomial at t
  # The polynomial.functions package stores coefficients in increasing degree order
  # So we can use polynomial::as.function() or evaluate manually
  result <- as.function(poly_k)(t)
  
  return(result)
}


#' Scale Timepoints to [-1, 1]
#'
#' Linearly transforms timepoints to the interval [-1, 1] required for
#' Legendre polynomial evaluation.
#'
#' @param timepoints Numeric vector of timepoints (e.g., c(-15, 0, 30, 60, 90, 120))
#'
#' @return Numeric vector scaled to [-1, 1]
#'
#' @examples
#' timepoints <- c(-15, 0, 30, 60, 90, 120)
#' t_std <- scale_to_legendre_domain(timepoints)
#' # Returns: c(-1.0, -0.7778, -0.3333, 0.1111, 0.5556, 1.0)
#'
#' @export
scale_to_legendre_domain <- function(timepoints) {
  t_min <- min(timepoints, na.rm = TRUE)
  t_max <- max(timepoints, na.rm = TRUE)
  
  if (t_max == t_min) {
    stop("Cannot scale: all timepoints are identical")
  }
  
  t_std <- (timepoints - t_min) / (t_max - t_min) * 2 - 1
  return(t_std)
}


' Add Legendre Polynomial Basis Columns
#'
#' Adds Legendre polynomial basis functions (order 0-4) to long format data.
#' Also creates group-specific basis columns for multiple-group SEM framework.
#'
#' IMPORTANT: Requires legendre_polynomials.R to be sourced first.
#'
#' @param dat_long Long format data from prepare_ogtt_long()
#'
#' @return Data with additional columns:
#'   - t_std: scaled timepoints in [-1, 1]
#'   - leg_poly_0, leg_poly_1, ..., leg_poly_4: Legendre basis functions
#'   - leg_poly_X_ctrl, leg_poly_X_treat: group-specific basis columns
#'   - model_timepoint: converted to factor
#'
#' @examples
#' # Source Legendre utilities first
#' source("R/legendre_polynomials.R")
#' baseline_long <- baseline_long |> add_legendre_basis()
#'
#' @export
add_legendre_basis <- function(dat_long) {
  # Check if legendre_polynomials.R has been sourced
  if (!exists("leg_poly", mode = "function")) {
    stop("Function leg_poly() not found. Please source 'R/legendre_polynomials.R' first.")
  }
  if (!exists("scale_to_legendre_domain", mode = "function")) {
    stop("Function scale_to_legendre_domain() not found. Please source 'R/legendre_polynomials.R' first.")
  }
  
  dat_long = dat_long %>% 
    dplyr::mutate(
      t_std = scale_to_legendre_domain(model_timepoint),
      leg_poly_0 = 1,
      leg_poly_1 = leg_poly(t_std, 1),
      leg_poly_2 = leg_poly(t_std, 2),
      leg_poly_3 = leg_poly(t_std, 3),
      leg_poly_4 = leg_poly(t_std, 4),
      # Group-specific columns for multiple-group SEM
      leg_poly_0_ctrl = leg_poly_0 * (1 - treat),
      leg_poly_1_ctrl = leg_poly_1 * (1 - treat),
      leg_poly_2_ctrl = leg_poly_2 * (1 - treat),
      leg_poly_3_ctrl = leg_poly_3 * (1 - treat),
      leg_poly_4_ctrl = leg_poly_4 * (1 - treat),
      leg_poly_0_treat = leg_poly_0 * treat,
      leg_poly_1_treat = leg_poly_1 * treat,
      leg_poly_2_treat = leg_poly_2 * treat,
      leg_poly_3_treat = leg_poly_3 * treat,
      leg_poly_4_treat = leg_poly_4 * treat,
      model_timepoint = factor(model_timepoint, levels = c("0", "30", "60", "90", "120"))
    )
  
  print("names dat_long:")
  print(names(dat_long)) #sanity check
  return(dat_long)
  
}

#' Create Legendre Polynomial Design Matrix
#'
#' Constructs the full design matrix (Lambda) for Legendre polynomial growth
#' curves up to degree max_degree.
#'
#' @param timepoints Numeric vector of timepoints
#' @param max_degree Integer, maximum polynomial degree (default 4)
#' @param include_intercept Logical, include intercept column (default TRUE)
#'
#' @return Matrix with dimensions length(timepoints) × (max_degree + 1)
#'   Columns: Intercept, P_1(t), P_2(t), ..., P_max_degree(t)
#'
#' @examples
#' timepoints <- c(-15, 0, 30, 60, 90, 120)
#' Lambda <- legendre_design_matrix(timepoints, max_degree = 4)
#' # Returns 6 × 5 matrix: [1, P_1(t), P_2(t), P_3(t), P_4(t)]
#'
#' @export
legendre_design_matrix <- function(timepoints, max_degree = 4, include_intercept = TRUE) {
  t_std <- scale_to_legendre_domain(timepoints)
  n <- length(timepoints)
  
  if (include_intercept) {
    ncol <- max_degree + 1
    Lambda <- matrix(NA, nrow = n, ncol = ncol)
    Lambda[, 1] <- 1  # Intercept
    
    for (k in 1:max_degree) {
      Lambda[, k + 1] <- leg_poly(t_std, k)
    }
    
    colnames(Lambda) <- c("Intercept", paste0("P", 1:max_degree))
  } else {
    ncol <- max_degree
    Lambda <- matrix(NA, nrow = n, ncol = ncol)
    
    for (k in 1:max_degree) {
      Lambda[, k] <- leg_poly(t_std, k)
    }
    
    colnames(Lambda) <- paste0("P", 1:max_degree)
  }
  
  rownames(Lambda) <- paste0("t", timepoints)
  return(Lambda)
}

