#' Simulate a Gaussian random field
#'
#' Generates a Gaussian random field (GRF) on a window `W` using an
#' exponential covariance function.
#'
#' @param W Observation window (an `owin` object).
#' @param params Named vector with:
#'   - `mu`: mean of the field
#'   - `var`: variance of the field
#'   - `scale`: correlation scale
#'
#' @return A pixel image (`im`) representing the simulated GRF.
#' @export
#'
#' @examples
#' library(spatstat.geom)
#'
#' W <- owin(c(0, 10), c(0, 10))
#' g <- sim_grf(W, c(mu = 0, var = 0.2, scale = 1))
#' plot(g)
sim_grf <- function(
    W = spatstat.geom::owin(),
    params = c(mu = 0, scale = 0.1, var = 0.1)
) {
  # Check required parameters are present
  if (!all(c("mu", "var", "scale") %in% names(params))) {
    stop("params must contain mu, var, and scale.", call. = FALSE)
  }

  spatstat.random::rGRFexpo(
    W = W,
    mu = params[["mu"]],
    var = params[["var"]],
    scale = params[["scale"]]
  )
}

#' Simulate an inhomogeneous Poisson process
#'
#' Simulates a Poisson point process with intensity given by a pixel image.
#'
#' @param lambda_im Intensity image (`im` object).
#' @param W Observation window.
#'
#' @return A point pattern (`ppp`).
#' @export
#'
#' @examples
#' library(spatstat.geom)
#'
#' W <- owin(c(0, 10), c(0, 10))
#' lambda <- as.im(function(x, y) 50, W)
#' X <- sim_pp(lambda, W)
#' plot(X)
sim_pp <- function(lambda_im, W) {
  spatstat.random::rpoispp(lambda_im, W = W)
}

#' Simulate a log-Gaussian Cox process
#'
#' This function:
#' 1. simulates a Gaussian random field,
#' 2. exponentiates it to get an intensity field,
#' 3. simulates a Poisson process from that intensity.
#'
#' You can supply either:
#' - (mu, sigma, s), or
#' - (lambda, phi, s)
#'
#' @param W Observation window.
#' @param mu Mean of the latent Gaussian random field.
#' @param sigma Standard deviation of the latent Gaussian random field.
#' @param s Correlation scale.
#' @param lambda Mean intensity.
#' @param phi Clustering strength, defined as g(0) - 1.
#'
#' @return A list with:
#' \describe{
#'   \item{pp}{Simulated point pattern.}
#'   \item{intensity}{Simulated intensity field.}
#'   \item{grf}{Latent Gaussian random field.}
#'   \item{params}{Parameters used in the simulation.}
#' }
sim_lgcp <- function(W = spatstat.geom::owin(),
                     mu = NULL, sigma = NULL, s = NULL,
                     lambda = NULL, phi = NULL) {
  # Convert the parameters to standard parameterization if necessary
  pars <- lgcp_params(
    mu = mu, sigma = sigma, s = s,
    lambda = lambda, phi = phi
  )

  grf <- sim_grf(
    params = c(mu = pars$mu, scale = pars$s, var = pars$var),
    W = W
  )

  # Exponentiate the GRF
  intensity <- exp(grf)

  # Generate the point pattern
  pp <- sim_pp(lambda_im = intensity, W = W)

  list(
    pp = pp,
    intensity = intensity,
    grf = grf,
    params = pars
  )
}

#' Simulate a CSCP
#'
#' This model has intensity
#'   Lambda(u) = mu + sum_{i=1}^df Z_i(u)^2
#'
#' where the Z_i are independent mean-zero Gaussian random fields.
#'
#' You can supply either:
#' - (mu, sigma, s), or
#' - (lambda, phi, s)
#'
#' @param W Observation window.
#' @param mu Baseline intensity.
#' @param sigma Standard deviation of each latent Gaussian random field.
#' @param s Correlation scale.
#' @param lambda Mean intensity.
#' @param phi Clustering strength, defined as g(0) - 1.
#' @param df Number of GRFs in the chi-square part.
#'
#' @return A list with:
#' \describe{
#'   \item{pp}{Simulated point pattern.}
#'   \item{intensity}{Simulated intensity field.}
#'   \item{grf}{List of latent Gaussian random fields.}
#'   \item{params}{Parameters used in the simulation.}
#' }
sim_cscp <- function(W = spatstat.geom::owin(),
                      mu = NULL, sigma = NULL, s = NULL,
                      lambda = NULL, phi = NULL,
                      df = 1) {
  # Convert the parameters to standard parameterization if necessary
  pars <- cscp_params(
    mu = mu, sigma = sigma, s = s,
    lambda = lambda, phi = phi,
    df = df
  )

  # Simulate one GRF for each component, defaults to single component
  grf <- replicate(
    pars$df,
    sim_grf(
      params = c(mu = 0, scale = pars$s, var = pars$var),
      W = W
    ),
    simplify = FALSE
  )

  # Square each component and sum
  intensity <- Reduce(`+`, lapply(grf, function(g) g^2)) + pars$mu

  # Generate the point pattern
  pp <- sim_pp(lambda_im = intensity, W = W)

  list(
    pp = pp,
    intensity = intensity,
    grf = grf,
    params = pars
  )
}

#' Simulate a non-shifted, non-central CSCP
#'
#' This model has intensity
#'   Lambda(u) = sum_{i=1}^df Z_i(u)^2
#'
#' where the Z_i are independent Gaussian random fields with mean `mu`.
#'
#' You can supply either:
#' - (mu, sigma, s), or
#' - (lambda, phi, s)
#'
#' @param W Observation window.
#' @param mu Mean of each latent Gaussian random field.
#' @param sigma Standard deviation of each latent Gaussian random field.
#' @param s Correlation scale.
#' @param lambda Mean intensity.
#' @param phi Clustering strength, defined as g(0) - 1.
#' @param df Number of GRFs in the chi-square part.
#'
#' @return A list with:
#' \describe{
#'   \item{pp}{Simulated point pattern.}
#'   \item{intensity}{Simulated intensity field.}
#'   \item{grf}{List of latent Gaussian random fields.}
#'   \item{params}{Parameters used in the simulation.}
#' }
#' @export
sim_cscp_nc <- function(W = spatstat.geom::owin(),
                        mu = NULL, sigma = NULL, s = NULL,
                        lambda = NULL, phi = NULL,
                        df = 1) {

  pars <- cscp_nc_params(
    mu = mu, sigma = sigma, s = s,
    lambda = lambda, phi = phi,
    df = df
  )

  grf <- replicate(
    pars$df,
    sim_grf(
      params = c(mu = pars$mu, scale = pars$s, var = pars$var),
      W = W
    ),
    simplify = FALSE
  )

  intensity <- Reduce(`+`, lapply(grf, function(g) g^2))

  pp <- sim_pp(lambda_im = intensity, W = W)

  list(
    pp = pp,
    intensity = intensity,
    grf = grf,
    params = pars
  )
}
