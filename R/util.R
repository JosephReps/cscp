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
