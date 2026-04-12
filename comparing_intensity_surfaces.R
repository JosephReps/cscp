# Settings


set.seed(3)

W <- spatstat.geom::owin(c(0, 1), c(0, 1))

nrow_panel <- 5
ncol_panel <- 5
nsim <- nrow_panel * ncol_panel

# Common model parameters
lambda <- 1
phi    <- 2

# Model-specific scales
s_lgcp    <- 0.15
s_cscp    <- 0.21
s_cscp_nc <- 0.21

# For CSCP models
df_cscp    <- 1
df_cscp_nc <- 1

# Plot on log scale?
log_scale <- FALSE

# PCF plot range
r_max <- 0.5
n_r   <- 500

# True PCF functions


pcf_lgcp <- function(r, lambda = NULL, phi = NULL, s = NULL,
                     mu = NULL, sigma = NULL) {
  pars <- lgcp_params(
    mu = mu, sigma = sigma, s = s,
    lambda = lambda, phi = phi
  )

  exp(pars$var * exp(-r / pars$s))
}

pcf_cscp <- function(r, lambda = NULL, phi = NULL, s = NULL,
                     mu = NULL, sigma = NULL, df = 1) {
  pars <- cscp_params(
    mu = mu, sigma = sigma, s = s,
    lambda = lambda, phi = phi,
    df = df
  )

  1 + pars$phi * exp(-2 * r / pars$s)
}

pcf_cscp_nc <- function(r, lambda = NULL, phi = NULL, s = NULL,
                        mu = NULL, sigma = NULL, df = 1) {
  pars <- cscp_nc_params(
    mu = mu, sigma = sigma, s = s,
    lambda = lambda, phi = phi,
    df = df
  )

  rho <- exp(-r / pars$s)

  1 + (4 * pars$mu^2 * pars$var * rho + 2 * pars$var^2 * rho^2) /
    (pars$df * (pars$mu^2 + pars$var)^2)
}

# Helper: transform intensity for plotting

transform_im <- function(im, log_scale = TRUE) {
  if (!log_scale) return(im)

  out <- im
  out$v <- log(out$v)
  out
}


# Helper: simulate many intensity fields

sim_intensity_list <- function(model = c("lgcp", "cscp", "cscp_nc"),
                               nsim,
                               W,
                               lambda,
                               phi,
                               s,
                               df = 1,
                               log_scale = TRUE) {
  model <- match.arg(model)

  out <- vector("list", nsim)

  for (i in seq_len(nsim)) {
    sim <- switch(
      model,
      lgcp = sim_lgcp(W = W, lambda = lambda, phi = phi, s = s),
      cscp = sim_cscp(W = W, lambda = lambda, phi = phi, s = s, df = df),
      cscp_nc = sim_cscp_nc(W = W, lambda = lambda, phi = phi, s = s, df = df)
    )

    out[[i]] <- transform_im(sim$intensity, log_scale = log_scale)
  }

  out
}


# Simulate fields

lgcp_list <- sim_intensity_list(
  model = "lgcp",
  nsim = nsim,
  W = W,
  lambda = lambda,
  phi = phi,
  s = s_lgcp,
  log_scale = log_scale
)

cscp_list <- sim_intensity_list(
  model = "cscp",
  nsim = nsim,
  W = W,
  lambda = lambda,
  phi = phi,
  s = s_cscp,
  df = df_cscp,
  log_scale = log_scale
)

cscp_nc_list <- sim_intensity_list(
  model = "cscp_nc",
  nsim = nsim,
  W = W,
  lambda = lambda,
  phi = phi,
  s = s_cscp_nc,
  df = df_cscp_nc,
  log_scale = log_scale
)


# Use common colour scale within each model

zlim_lgcp <- range(unlist(lapply(lgcp_list, function(x) x$v)), finite = TRUE)
zlim_cscp <- range(unlist(lapply(cscp_list, function(x) x$v)), finite = TRUE)
zlim_nc   <- range(unlist(lapply(cscp_nc_list, function(x) x$v)), finite = TRUE)


# UNCOMMENT BELOW FOR SHARED COLOUR SCALES
zlim_all <- range(
  c(
    unlist(lapply(lgcp_list, function(x) x$v)),
    unlist(lapply(cscp_list, function(x) x$v)),
    unlist(lapply(cscp_nc_list, function(x) x$v))
  ),
  finite = TRUE
)

zlim_lgcp <- zlim_all
zlim_cscp <- zlim_all
zlim_nc   <- zlim_all


# Figure 1: three 5x5 panels

old_par <- par(no.readonly = TRUE)

layout(matrix(1:3, nrow = 1), widths = c(1, 1, 1))

draw_panel <- function(im_list, zlim, nrow_panel, ncol_panel) {
  par(mfrow = c(nrow_panel, ncol_panel),
      mar = c(0, 0, 0, 0),
      oma = c(0, 0, 0, 0),
      xaxs = "i",
      yaxs = "i")

  for (im in im_list) {
    plot(
      im,
      ribbon = FALSE,
      axes = FALSE,
      box = FALSE,
      main = "",
      useRaster = TRUE,
      zlim = zlim
    )
  }
}

draw_panel(lgcp_list, zlim_lgcp, nrow_panel, ncol_panel)
draw_panel(cscp_list, zlim_cscp, nrow_panel, ncol_panel)
draw_panel(cscp_nc_list, zlim_nc, nrow_panel, ncol_panel)

par(old_par)


# Figure 2: theoretical PCF curves


r <- seq(0, r_max, length.out = n_r)

g_lgcp <- pcf_lgcp(
  r = r,
  lambda = lambda,
  phi = phi,
  s = s_lgcp
)

g_cscp <- pcf_cscp(
  r = r,
  lambda = lambda,
  phi = phi,
  s = s_cscp,
  df = df_cscp
)

g_cscp_nc <- pcf_cscp_nc(
  r = r,
  lambda = lambda,
  phi = phi,
  s = s_cscp_nc,
  df = df_cscp_nc
)

par(mar = c(3, 3, 1, 1))

plot(
  r, g_lgcp,
  type = "l",
  lwd = 2,
  xlab = "r",
  ylab = "g(r)",
  ylim = range(c(g_lgcp, g_cscp, g_cscp_nc))
)

lines(r, g_cscp, lwd = 2, lty = 2)
lines(r, g_cscp_nc, lwd = 2, lty = 3)

legend(
  "topright",
  legend = c("LGCP", "CSCP", "CSCP_NC"),
  lwd = 2,
  lty = c(1, 2, 3),
  bty = "n"
)
