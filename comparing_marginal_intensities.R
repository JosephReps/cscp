set.seed(2)

W <- spatstat.geom::owin(c(0, 1), c(0, 1))

lgcp    <- sim_lgcp(W, lambda = 1, phi = 1, s = 0.07)
cscp_s  <- sim_cscp(W, lambda = 1, phi = 1, s = 0.05, df = 1)      # shifted
cscp_nc <- sim_cscp_nc(W, lambda = 1, phi = 1, s = 0.05, df = 1)   # non-shifted non-central

z_lgcp    <- log(as.vector(lgcp$intensity$v))
z_cscp_s  <- log(as.vector(cscp_s$intensity$v))
z_cscp_nc <- log(as.vector(cscp_nc$intensity$v))

# z_lgcp    <- exp(z_lgcp)
# z_cscp_s  <- exp(z_cscp_s)
# z_cscp_nc <- exp(z_cscp_nc)

xlim <- range(c(z_lgcp, z_cscp_s, z_cscp_nc), finite = TRUE)

## Histogram overlayed
all_vals <- c(z_lgcp, z_cscp_s, z_cscp_nc)

# breaks <- seq(min(all_vals), max(all_vals), length.out = 60)
#
# hist(z_lgcp, breaks = breaks, col = rgb(0,0,0,0.4), main = "", xlab = "")
# hist(z_cscp_s, breaks = breaks, col = rgb(1,0,0,0.4), add = TRUE)
# hist(z_cscp_nc, breaks = breaks, col = rgb(0,0,1,0.4), add = TRUE)

## Density plot

# plot(density(z_lgcp), col = "black", lwd = 2, main = "", xlab = "", ylim = c(0,1.2))
# lines(density(z_cscp_s), col = "red", lwd = 2)
# lines(density(z_cscp_nc), col = "blue", lwd = 2)
# legend("topright",
#        legend = c("LGCP", "CSCP (shifted)", "CSCP (non-central)"),
#        col = c("black", "red", "blue"),
#        lwd = 2,
#        bty = "n")


################## Theory + empirical marginal intensity density

p_lgcp   <- lgcp$params
p_cscp_s <- cscp_s$params
p_cscp_nc <- cscp_nc$params

# Theoretical densities (log scale)

d_lgcp_log <- function(x, pars) {
  dnorm(x, mean = pars$mu, sd = pars$sigma)
}

d_cscp_shifted_log <- function(x, pars) {
  mu <- pars$mu
  sigma <- pars$sigma

  out <- numeric(length(x))

  valid <- x > log(mu)
  z <- exp(x[valid]) - mu

  out[valid] <- (1 / (sigma * sqrt(2*pi))) *
    (1 / sqrt(z)) *
    exp(-z / (2 * sigma^2)) *
    exp(x[valid])

  out
}

d_cscp_nc_log <- function(x, pars) {
  mu <- pars$mu
  sigma <- pars$sigma

  z <- exp(x)
  root <- sqrt(z)

  dens <- (1 / (sigma * sqrt(2*pi*z))) * (
    exp(-(root - mu)^2 / (2 * sigma^2)) +
      exp(-(root + mu)^2 / (2 * sigma^2))
  )

  dens * z * 0.5
}

# grid
x_grid <- seq(min(all_vals), max(all_vals), length.out = 500)

# empirical
plot(density(z_lgcp), col = "black", lwd = 2, lty = 2,
     main = "", xlab = "", ylim = c(0, 1))
lines(density(z_cscp_s), col = "red", lwd = 2, lty = 2)
lines(density(z_cscp_nc), col = "blue", lwd = 2, lty = 2)

# theoretical
lines(x_grid, d_lgcp_log(x_grid, p_lgcp), col = "black", lwd = 2)
lines(x_grid, d_cscp_shifted_log(x_grid, p_cscp_s), col = "red", lwd = 2)
lines(x_grid, d_cscp_nc_log(x_grid, p_cscp_nc), col = "blue", lwd = 2)

legend("topright",
       legend = c("LGCP (emp)", "CSCP shifted (emp)", "CSCP nc (emp)",
                  "LGCP (theory)", "CSCP shifted (theory)", "CSCP nc (theory)"),
       col = c("black","red","blue","black","red","blue"),
       lty = c(2,2,2,1,1,1),
       lwd = 2,
       bty = "n")







