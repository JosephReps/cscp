#' Convert LGCP parameters
#'
#' You can supply either:
#' - (mu, sigma, s), or
#' - (lambda, phi, s)
#'
#' For the LGCP:
#' - lambda = exp(mu + sigma^2 / 2)
#' - phi    = exp(sigma^2) - 1
#'
#' @param mu Mean of the latent Gaussian random field.
#' @param sigma Standard deviation of the latent Gaussian random field.
#' @param s Correlation scale.
#' @param lambda Mean intensity.
#' @param phi Clustering strength, defined as g(0) - 1.
#'
#' @return A named list with mu, sigma, var, lambda, phi, s.
lgcp_params <- function(mu = NULL, sigma = NULL, s = NULL,
                        lambda = NULL, phi = NULL) {

  # Make sure we are providing exactly one valid parameterization
  using_standard <- !is.null(mu) || !is.null(sigma)
  using_alt <- !is.null(lambda) || !is.null(phi)

  if (using_standard && using_alt) {
    stop("Supply either (mu, sigma, s) or (lambda, phi, s), not both.", call. = FALSE)
  }

  if (!using_standard && !using_alt) {
    stop("Supply either (mu, sigma, s) or (lambda, phi, s).", call. = FALSE)
  }

  # If standard parameterization provided
  if (using_standard) {
    # Make sure everything is valid
    if (is.null(mu) || is.null(sigma) || is.null(s)) {
      stop("For LGCP, (mu, sigma, s) must all be supplied.", call. = FALSE)
    }
    stopifnot(is.numeric(mu), length(mu) == 1, is.finite(mu))
    stopifnot(is.numeric(sigma), length(sigma) == 1, is.finite(sigma), sigma >= 0)
    stopifnot(is.numeric(s), length(s) == 1, is.finite(s), s > 0)

    # Then calculate the cluster-strength parameters
    var <- sigma^2
    lambda <- exp(mu + var / 2)
    phi <- exp(var) - 1
  }

  # If using the clustering-strength parameterization
  if (using_alt) {
    # Make sure everything is valid
    if (is.null(lambda) || is.null(phi) || is.null(s)) {
      stop("For LGCP, (lambda, phi, s) must all be supplied.", call. = FALSE)
    }
    stopifnot(is.numeric(lambda), length(lambda) == 1, is.finite(lambda), lambda > 0)
    stopifnot(is.numeric(phi), length(phi) == 1, is.finite(phi), phi >= 0)
    stopifnot(is.numeric(s), length(s) == 1, is.finite(s), s > 0)

    # Calculate the standard parameters
    var <- log(1 + phi)
    sigma <- sqrt(var)
    mu <- log(lambda) - var / 2
  }

  # Return them all, let user decide what to do
  list(
    mu = unname(mu),
    sigma = unname(sigma),
    var = unname(var),
    lambda = unname(lambda),
    phi = unname(phi),
    s = unname(s)
  )
}

#' Convert CSCP parameters
#'
#' CSCP has intensity
#'   Lambda(u) = mu + sum_{i=1}^df Z_i(u)^2
#'
#' where each Z_i is a mean-zero Gaussian random field with standard
#' deviation sigma.
#'
#' You can supply either:
#' - (mu, sigma, s), or
#' - (lambda, phi, s)
#'
#' For this model:
#' - lambda = mu + df * sigma^2
#' - phi    = 2 * df * sigma^4 / lambda^2
#'
#' @param mu Baseline intensity.
#' @param sigma Standard deviation of each latent Gaussian random field.
#' @param s Correlation scale.
#' @param lambda Mean intensity.
#' @param phi Clustering strength, defined as g(0) - 1.
#' @param df Number of GRFs in the chi-square part.
#'
#' @return A named list with mu, sigma, var, lambda, phi, s, df.
cscp_params <- function(mu = NULL, sigma = NULL, s = NULL,
                         lambda = NULL, phi = NULL,
                         df = 1) {
  # CSCP slightly mroe tricky than LGCP, as we have to account for the
  # possibility of multiple components, but still straightforward

  stopifnot(is.numeric(df), length(df) == 1, is.finite(df), df >= 1)
  df <- as.integer(df)

  # Make sure we are providing exactly one valid parameterization
  using_standard <- !is.null(mu) || !is.null(sigma)
  using_alt <- !is.null(lambda) || !is.null(phi)

  if (using_standard && using_alt) {
    stop("Supply either (mu, sigma, s) or (lambda, phi, s), not both.", call. = FALSE)
  }

  if (!using_standard && !using_alt) {
    stop("Supply either (mu, sigma, s) or (lambda, phi, s).", call. = FALSE)
  }

  # If standard parameterization provided
  if (using_standard) {
    # Make sure everything is valid
    if (is.null(mu) || is.null(sigma) || is.null(s)) {
      stop("For CSCP, (mu, sigma, s) must all be supplied.", call. = FALSE)
    }
    stopifnot(is.numeric(mu), length(mu) == 1, is.finite(mu), mu >= 0)
    stopifnot(is.numeric(sigma), length(sigma) == 1, is.finite(sigma), sigma >= 0)
    stopifnot(is.numeric(s), length(s) == 1, is.finite(s), s > 0)

    # Calculate the cluster-strength parameters
    var <- sigma^2
    lambda <- mu + df * var
    phi <- if (lambda == 0) 0 else 2 * df * var^2 / lambda^2
  }

  # If using the clustering-strength parameterization
  if (using_alt) {
    # Make sure everything is valid
    if (is.null(lambda) || is.null(phi) || is.null(s)) {
      stop("For CSCP, (lambda, phi, s) must all be supplied.", call. = FALSE)
    }
    stopifnot(is.numeric(lambda), length(lambda) == 1, is.finite(lambda), lambda > 0)
    stopifnot(is.numeric(phi), length(phi) == 1, is.finite(phi), phi >= 0)
    stopifnot(is.numeric(s), length(s) == 1, is.finite(s), s > 0)

    if (phi > 2 / df) {
      # stop("For CSCP, phi must be <= 2 / df.", call. = FALSE)
      warning("For CSCP, phi must be <= 2 / df.")
    }

    # Calculate the standard parameters
    var <- lambda * sqrt(phi / (2 * df))
    sigma <- sqrt(var)
    mu <- lambda - df * var
  }

  # Return them all
  list(
    mu = unname(mu),
    sigma = unname(sigma),
    var = unname(var),
    lambda = unname(lambda),
    phi = unname(phi),
    s = unname(s),
    df = df
  )
}

#' Convert non-shifted, non-central CSCP parameters
#'
#' This model has intensity
#'   Lambda(u) = sum_{i=1}^df Z_i(u)^2
#'
#' where each Z_i is a Gaussian random field with mean `mu` and standard
#' deviation `sigma`.
#'
#' You can supply either:
#' - (mu, sigma, s), or
#' - (lambda, phi, s)
#'
#' For this model:
#' - lambda = df * (mu^2 + sigma^2)
#' - phi    = [4 mu^2 sigma^2 + 2 sigma^4] / [df (mu^2 + sigma^2)^2]
#'
#' @param mu Mean of each latent Gaussian random field.
#' @param sigma Standard deviation of each latent Gaussian random field.
#' @param s Correlation scale.
#' @param lambda Mean intensity.
#' @param phi Clustering strength, defined as g(0) - 1.
#' @param df Number of GRFs in the chi-square part.
#'
#' @return A named list with mu, sigma, var, lambda, phi, s, df.
#' @export
cscp_nc_params <- function(mu = NULL, sigma = NULL, s = NULL,
                           lambda = NULL, phi = NULL,
                           df = 1) {

  stopifnot(is.numeric(df), length(df) == 1, is.finite(df), df >= 1)
  df <- as.integer(df)

  using_standard <- !is.null(mu) || !is.null(sigma)
  using_alt <- !is.null(lambda) || !is.null(phi)

  if (using_standard && using_alt) {
    stop("Supply either (mu, sigma, s) or (lambda, phi, s), not both.", call. = FALSE)
  }

  if (!using_standard && !using_alt) {
    stop("Supply either (mu, sigma, s) or (lambda, phi, s).", call. = FALSE)
  }

  # Standard parameterization: (mu, sigma, s)
  if (using_standard) {
    if (is.null(mu) || is.null(sigma) || is.null(s)) {
      stop("For non-central CSCP, (mu, sigma, s) must all be supplied.", call. = FALSE)
    }

    stopifnot(is.numeric(mu), length(mu) == 1, is.finite(mu), mu >= 0)
    stopifnot(is.numeric(sigma), length(sigma) == 1, is.finite(sigma), sigma >= 0)
    stopifnot(is.numeric(s), length(s) == 1, is.finite(s), s > 0)

    var <- sigma^2
    A <- mu^2 + var

    lambda <- df * A
    phi <- if (A == 0) 0 else (4 * mu^2 * var + 2 * var^2) / (df * A^2)
  }

  # Alternative parameterization: (lambda, phi, s)
  if (using_alt) {
    if (is.null(lambda) || is.null(phi) || is.null(s)) {
      stop("For non-central CSCP, (lambda, phi, s) must all be supplied.", call. = FALSE)
    }

    stopifnot(is.numeric(lambda), length(lambda) == 1, is.finite(lambda), lambda > 0)
    stopifnot(is.numeric(phi), length(phi) == 1, is.finite(phi), phi >= 0)
    stopifnot(is.numeric(s), length(s) == 1, is.finite(s), s > 0)

    if (phi > 2 / df) {
      stop("For non-central CSCP, phi must be <= 2 / df.", call. = FALSE)
    }

    A <- lambda / df
    q <- sqrt(1 - df * phi / 2)

    var <- A * (1 - q)
    sigma <- sqrt(var)

    mu2 <- A * q
    mu <- sqrt(mu2)
  }

  list(
    mu = unname(mu),
    sigma = unname(sigma),
    var = unname(var),
    lambda = unname(lambda),
    phi = unname(phi),
    s = unname(s),
    df = df
  )
}

#' Match an LGCP scale to a CSCP scale
#'
#' Numerically finds the LGCP correlation scale that best matches the
#' theoretical CSCP pair correlation function, for fixed clustering strength
#' `phi`, by minimizing the sup-norm discrepancy over a grid in
#' \eqn{x \in [0,1]}.
#'
#' The matching is based on the theoretical pair correlation functions
#'
#' \deqn{
#' g_{CSCP}(r) = 1 + \phi \exp(-2r / s_C)
#' }
#'
#' and
#'
#' \deqn{
#' g_{LGCP}(r) = (1+\phi)^{\exp(-r / s_L)}.
#' }
#'
#' Writing \eqn{\alpha = s_L / s_C}, the function numerically minimizes
#'
#' \deqn{
#' \sup_{x \in [0,1]}
#' \left| \exp\{ \log(1+\phi)\, x^{1/(2\alpha)} \} - 1 - \phi x \right|.
#' }
#'
#' @param phi Numeric scalar or vector of clustering strengths.
#'   Must satisfy `phi >= 0`.
#' @param s_cscp Numeric scalar or vector of CSCP correlation scales.
#'   Must be strictly positive.
#' @param x_grid Numeric grid on `[0, 1]` used to approximate the sup norm.
#' @param alpha_lower Lower bound for the scale ratio
#'   \eqn{\alpha = s_L / s_C}.
#' @param alpha_upper Upper bound for the scale ratio
#'   \eqn{\alpha = s_L / s_C}.
#' @param return_full Logical; if `FALSE`, returns the matched LGCP scale(s).
#'   If `TRUE`, returns a data frame with matching details.
#'
#' @return If `return_full = FALSE`, a numeric vector of matched LGCP scales.
#'
#' If `return_full = TRUE`, a data frame with columns:
#' \describe{
#'   \item{phi}{Clustering strength.}
#'   \item{s_cscp}{Input CSCP scale.}
#'   \item{alpha_opt}{Optimal scale ratio `s_lgcp / s_cscp`.}
#'   \item{s_lgcp}{Matched LGCP scale.}
#'   \item{p_opt}{Equivalent exponent `1 / (2 * alpha_opt)`.}
#'   \item{match_error}{Approximated minimum sup-norm discrepancy.}
#' }
#'
#' @export
#'
#' @examples
#' match_lgcp_scale(phi = 1, s_cscp = 0.05)
#'
#' match_lgcp_scale(
#'   phi = c(0.1, 0.8, 1.5),
#'   s_cscp = 0.05,
#'   return_full = TRUE
#' )
match_lgcp_scale <- function(phi,
                             s_cscp,
                             x_grid = seq(0, 1, length.out = 4000),
                             alpha_lower = 0.05,
                             alpha_upper = 5,
                             return_full = FALSE) {

  if (!is.numeric(phi) || any(!is.finite(phi)) || any(phi < 0)) {
    stop("phi must be a numeric vector of finite non-negative values.",
         call. = FALSE)
  }

  if (!is.numeric(s_cscp) || any(!is.finite(s_cscp)) || any(s_cscp <= 0)) {
    stop("s_cscp must be a numeric vector of finite positive values.",
         call. = FALSE)
  }

  if (!is.numeric(x_grid) || any(!is.finite(x_grid)) ||
      any(x_grid < 0) || any(x_grid > 1)) {
    stop("x_grid must be a numeric vector with values in [0, 1].",
         call. = FALSE)
  }

  if (!is.numeric(alpha_lower) || length(alpha_lower) != 1 ||
      !is.finite(alpha_lower) || alpha_lower <= 0) {
    stop("alpha_lower must be a single finite positive number.",
         call. = FALSE)
  }

  if (!is.numeric(alpha_upper) || length(alpha_upper) != 1 ||
      !is.finite(alpha_upper) || alpha_upper <= alpha_lower) {
    stop("alpha_upper must be a single finite number greater than alpha_lower.",
         call. = FALSE)
  }

  n <- max(length(phi), length(s_cscp))

  if (!(length(phi) %in% c(1, n) && length(s_cscp) %in% c(1, n))) {
    stop("phi and s_cscp must have compatible lengths.",
         call. = FALSE)
  }

  phi <- rep_len(phi, n)
  s_cscp <- rep_len(s_cscp, n)

  match_one <- function(phi_i, s_i) {
    obj_alpha <- function(alpha) {
      p <- 1 / (2 * alpha)
      cval <- log(1 + phi_i)
      max(abs(exp(cval * x_grid^p) - 1 - phi_i * x_grid))
    }

    opt <- stats::optimize(
      f = obj_alpha,
      interval = c(alpha_lower, alpha_upper)
    )

    alpha_opt <- opt$minimum
    s_lgcp <- alpha_opt * s_i

    data.frame(
      phi = phi_i,
      s_cscp = s_i,
      alpha_opt = alpha_opt,
      s_lgcp = s_lgcp,
      p_opt = 1 / (2 * alpha_opt),
      match_error = opt$objective
    )
  }

  out <- do.call(
    rbind,
    Map(match_one, phi_i = phi, s_i = s_cscp)
  )

  rownames(out) <- NULL

  if (return_full) {
    return(out)
  }

  out$s_lgcp
}
