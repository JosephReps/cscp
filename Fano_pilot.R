# ============================================================
# Pilot test: multi-scale quadrat count variability
# LGCP vs CSCP
# ============================================================

rm(list = ls())

library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)

set.seed(1)

# ============================================================
# Window
# ============================================================

W <- owin(c(0, 1), c(0, 1))

# ============================================================
# Parameters
# ============================================================

lambda_target <- 300
phi_target    <- 1.5
s_target      <- 0.05
s_target_LGCP      <- s_target * 0.68

n_rep <- 100

# Quadrat side lengths to try
cell_sizes <- seq(0.02, 0.4, by = 0.02)

# ============================================================
# Helper functions
# ============================================================

# LGCP parameter mapping:
# g(0) - 1 = exp(sigma2) - 1 = phi  => sigma2 = log(1 + phi)
# lambda = exp(mu + sigma2/2)       => mu = log(lambda) - sigma2/2
sim_lgcp_matched <- function(W, lambda, phi, s) {
  sigma2 <- log(1 + phi)
  mu     <- log(lambda) - sigma2 / 2

  X <- rLGCP(
    model = "exponential",
    mu    = mu,
    var   = sigma2,
    scale = s,
    win   = W
  )

  X
}

# CSCP simulation via discretised Gaussian field:
# Lambda(u) = mu + Z(u)^2
# with Var(Z(u)) = sigma2 and
# phi = 2 sigma2^2 / lambda^2  => sigma2 = lambda * sqrt(phi / 2)
# mu = lambda - sigma2
#
# This requires phi <= 2 if mu >= 0.
sim_cscp_matched <- function(W, lambda, phi, s,
                             nx = 128, ny = 128) {

  sigma2 <- lambda * sqrt(phi / 2)
  mu     <- lambda - sigma2

  if (mu < 0) {
    stop("CSCP parameterisation gives mu < 0. Need phi <= 2.")
  }

  # simulate Gaussian random field on grid
  xseq <- seq(W$xrange[1], W$xrange[2], length.out = nx)
  yseq <- seq(W$yrange[1], W$yrange[2], length.out = ny)

  # random field simulation via RandomFields-style approximation is not in spatstat,
  # so here we use rLGCP's latent Gaussian field trick indirectly:
  # simulate GRF Y with mean 0 by simulating LGCP latent field then undoing exp is not available,
  # so instead use a correlated Gaussian image from blur of white noise as a quick pilot.
  #
  # This is only a pilot approximation, not your final thesis-grade simulator.

  Zmat <- matrix(rnorm(nx * ny, mean = 0, sd = sqrt(sigma2)), nrow = nx, ncol = ny)

  # crude spatial smoothing to induce correlation
  Zim <- im(t(Zmat), xcol = xseq, yrow = yseq)
  Zim <- blur(Zim, sigma = s / 2)

  # re-standardise after smoothing
  zvals <- as.vector(Zim$v)
  zvals <- (zvals - mean(zvals, na.rm = TRUE)) / sd(zvals, na.rm = TRUE)
  Zim$v <- matrix(zvals * sqrt(sigma2), nrow = nrow(Zim$v), ncol = ncol(Zim$v))

  Lambda <- eval.im(mu + Zim^2)

  X <- rpoispp(Lambda)

  X
}

# Compute quadrat summary statistics at one cell size
quadrat_stats_one_scale <- function(X, cell_size) {
  nx <- max(1, floor(diff(X$window$xrange) / cell_size))
  ny <- max(1, floor(diff(X$window$yrange) / cell_size))

  q <- quadratcount(X, nx = nx, ny = ny)
  counts <- as.vector(q)

  m <- mean(counts)
  v <- var(counts)

  data.frame(
    cell_size = cell_size,
    nx = nx,
    ny = ny,
    mean_count = m,
    var_count = v,
    fano = ifelse(m > 0, v / m, NA_real_),
    cv = ifelse(m > 0, sqrt(v) / m, NA_real_)
  )
}

# Compute stats across all scales for one pattern
quadrat_stats_all_scales <- function(X, cell_sizes) {
  bind_rows(lapply(cell_sizes, function(cs) quadrat_stats_one_scale(X, cs)))
}

# ============================================================
# Run pilot
# ============================================================

results <- vector("list", length = 2 * n_rep)
counter <- 1

for (i in seq_len(n_rep)) {
  cat("Simulation", i, "of", n_rep, "\n")

  X_lgcp <- sim_lgcp_matched(W, lambda_target, phi_target, s_target_LGCP)
  X_cscp <- sim_cscp_matched(W, lambda_target, phi_target, s_target)

  stats_lgcp <- quadrat_stats_all_scales(X_lgcp, cell_sizes) %>%
    mutate(model = "LGCP", rep = i)

  stats_cscp <- quadrat_stats_all_scales(X_cscp, cell_sizes) %>%
    mutate(model = "CSCP", rep = i)

  results[[counter]] <- stats_lgcp
  counter <- counter + 1

  results[[counter]] <- stats_cscp
  counter <- counter + 1
}

results_df <- bind_rows(results)

# ============================================================
# Aggregate
# ============================================================

summary_df <- results_df %>%
  group_by(model, cell_size) %>%
  summarise(
    mean_fano = mean(fano, na.rm = TRUE),
    sd_fano   = sd(fano, na.rm = TRUE),
    mean_cv   = mean(cv, na.rm = TRUE),
    sd_cv     = sd(cv, na.rm = TRUE),
    mean_var  = mean(var_count, na.rm = TRUE),
    mean_mean = mean(mean_count, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_df)

# ============================================================
# Plots
# ============================================================

p1 <- ggplot(summary_df, aes(x = cell_size, y = mean_fano, colour = model)) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    title = "Fano factor across quadrat scales",
    x = "Quadrat side length",
    y = "Mean Fano factor"
  ) +
  theme_minimal()

p2 <- ggplot(summary_df, aes(x = cell_size, y = mean_cv, colour = model)) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    title = "Coefficient of variation across quadrat scales",
    x = "Quadrat side length",
    y = "Mean CV"
  ) +
  theme_minimal()

p3 <- ggplot(summary_df, aes(x = cell_size, y = mean_var, colour = model)) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    title = "Variance of quadrat counts across scales",
    x = "Quadrat side length",
    y = "Mean variance"
  ) +
  theme_minimal()

print(p1)
print(p2)
print(p3)

# ============================================================
# Optional: inspect a single realisation
# ============================================================

X_lgcp_demo <- sim_lgcp_matched(W, lambda_target, phi_target, s_target_LGCP)
X_cscp_demo <- sim_cscp_matched(W, lambda_target, phi_target, s_target)

par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
plot(X_cscp_demo, main = "CSCP example")
plot(X_lgcp_demo, main = "LGCP example")
