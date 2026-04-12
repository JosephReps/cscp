#' Theoretical pair correlation function for an LGCP
#'
#' Computes the theoretical pair correlation function
#'
#' \deqn{
#' g(r) = \exp\left(\sigma^2 \exp(-r / s)\right)
#' }
#'
#' for a log-Gaussian Cox process with exponential covariance.
#'
#' You can supply either:
#' - (mu, sigma, s), or
#' - (lambda, phi, s)
#'
#' where
#' - \eqn{\lambda = \exp(\mu + \sigma^2 / 2)}
#' - \eqn{\phi = g(0) - 1 = \exp(\sigma^2) - 1}
#'
#' @param r Numeric vector of distances.
#' @param mu Mean of the latent Gaussian random field.
#' @param sigma Standard deviation of the latent Gaussian random field.
#' @param s Correlation scale.
#' @param lambda Mean intensity.
#' @param phi Clustering strength, defined as \eqn{g(0) - 1}.
#'
#' @return A numeric vector of the same length as `r`.
#' @export
#'
#' @examples
#' r <- seq(0, 1, length.out = 100)
#' g <- pcf_lgcp_theoretical(r, lambda = 1, phi = 1, s = 0.1)
#' plot(r, g, type = "l")
pcf_lgcp_theoretical <- function(r,
                                 mu = NULL, sigma = NULL, s = NULL,
                                 lambda = NULL, phi = NULL) {

  if (!is.numeric(r) || any(!is.finite(r)) || any(r < 0)) {
    stop("r must be a numeric vector of finite non-negative values.", call. = FALSE)
  }

  pars <- lgcp_params(
    mu = mu, sigma = sigma, s = s,
    lambda = lambda, phi = phi
  )

  exp(pars$var * exp(-r / pars$s))
}

#' Theoretical pair correlation function for a shifted CSCP
#'
#' Computes the theoretical pair correlation function
#'
#' \deqn{
#' g(r) = 1 + \phi \exp(-2r / s)
#' }
#'
#' for the shifted chi-square Cox process
#'
#' \deqn{
#' \Lambda(u) = \mu + \sum_{i=1}^{df} Z_i(u)^2,
#' }
#'
#' where the \eqn{Z_i} are independent mean-zero Gaussian random fields
#' with exponential correlation.
#'
#' You can supply either:
#' - (mu, sigma, s), or
#' - (lambda, phi, s)
#'
#' where
#' - \eqn{\lambda = \mu + df \sigma^2}
#' - \eqn{\phi = g(0) - 1 = 2 df \sigma^4 / \lambda^2}
#'
#' @param r Numeric vector of distances.
#' @param mu Baseline intensity.
#' @param sigma Standard deviation of each latent Gaussian random field.
#' @param s Correlation scale.
#' @param lambda Mean intensity.
#' @param phi Clustering strength, defined as \eqn{g(0) - 1}.
#' @param df Number of Gaussian random field components.
#'
#' @return A numeric vector of the same length as `r`.
#' @export
#'
#' @examples
#' r <- seq(0, 1, length.out = 100)
#' g <- pcf_cscp_theoretical(r, lambda = 1, phi = 0.5, s = 0.1, df = 1)
#' plot(r, g, type = "l")
pcf_cscp_theoretical <- function(r,
                                 mu = NULL, sigma = NULL, s = NULL,
                                 lambda = NULL, phi = NULL,
                                 df = 1) {

  if (!is.numeric(r) || any(!is.finite(r)) || any(r < 0)) {
    stop("r must be a numeric vector of finite non-negative values.", call. = FALSE)
  }

  pars <- cscp_params(
    mu = mu, sigma = sigma, s = s,
    lambda = lambda, phi = phi,
    df = df
  )

  1 + pars$phi * exp(-2 * r / pars$s)
}
