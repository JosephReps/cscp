# Pilot test: K-function comparison
rm(list = ls())
library(spatstat)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)

set.seed(1)
W <- owin(c(0, 1), c(0, 1))

# Parameters
lambda_target <- 300
phi_target    <- 1.5
s_target      <- 0.05
s_target_LGCP      <- 0.05 * 0.68

n_rep <- 100

# Grid of r values for K-function comparison
r_vals <- seq(0, 0.25, by = 0.005)

# Helper functionsssss
# LGCP parameter mapping:
# g(0) - 1 = exp(sigma2) - 1 = phi
# sigma2 = log(1 + phi)
# lambda = exp(mu + sigma2 / 2)
# mu = log(lambda) - sigma2 / 2
sim_lgcp_matched <- function(W, lambda, phi, s) {
  sigma2 <- log(1 + phi)
  mu     <- log(lambda) - sigma2 / 2

  rLGCP(
    model = "exponential",
    mu    = mu,
    var   = sigma2,
    scale = s,
    win   = W
  )
}

# Quick pilot CSCP approximation:
# Lambda(u) = mu + Z(u)^2
# phi = 2 sigma2^2 / lambda^2  => sigma2 = lambda * sqrt(phi / 2)
# mu = lambda - sigma2
sim_cscp_matched <- function(W, lambda, phi, s, nx = 128, ny = 128) {
  sigma2 <- lambda * sqrt(phi / 2)
  mu     <- lambda - sigma2

  if (mu < 0) {
    stop("CSCP parameterisation gives mu < 0. Need phi <= 2.")
  }

  xseq <- seq(W$xrange[1], W$xrange[2], length.out = nx)
  yseq <- seq(W$yrange[1], W$yrange[2], length.out = ny)

  Zmat <- matrix(rnorm(nx * ny, mean = 0, sd = sqrt(sigma2)), nrow = nx, ncol = ny)
  Zim  <- im(t(Zmat), xcol = xseq, yrow = yseq)

  # crude smoothing to induce correlation
  Zim <- blur(Zim, sigma = s / 2)

  # standardise back to target variance
  zvals <- as.vector(Zim$v)
  zvals <- (zvals - mean(zvals, na.rm = TRUE)) / sd(zvals, na.rm = TRUE)
  Zim$v <- matrix(
    zvals * sqrt(sigma2),
    nrow = nrow(Zim$v),
    ncol = ncol(Zim$v)
  )

  Lambda <- eval.im(mu + Zim^2)

  rpoispp(Lambda)
}

# Theoretical PCFs
pcf_lgcp_theoretical <- function(r, phi, s) {
  (1 + phi)^(exp(-r / s))
}

pcf_cscp_theoretical <- function(r, phi, s) {
  1 + phi * exp(-2 * r / s)
}

# Numerical K-function from PCF:
# K(r) = 2 pi \int_0^r t g(t) dt
k_from_pcf <- function(r_vals, gfun) {
  out <- numeric(length(r_vals))

  for (i in seq_along(r_vals)) {
    r <- r_vals[i]
    if (r == 0) {
      out[i] <- 0
    } else {
      integrand <- function(t) t * gfun(t)
      out[i] <- 2 * pi * integrate(integrand, lower = 0, upper = r,
                                   subdivisions = 1000L,
                                   rel.tol = 1e-6)$value
    }
  }

  out
}

# Extract empirical K and L-r on common r-grid
estimate_k_curve <- function(X, r_vals,
                             correction = "isotropic") {
  Kest_obj <- Kest(X, r = r_vals, correction = correction)

  # spatstat returns columns named like "iso", "trans", etc.
  col_name <- switch(
    correction,
    isotropic = "iso",
    translation = "trans",
    border = "border",
    "iso"
  )

  Khat <- Kest_obj[[col_name]]
  Lhat <- sqrt(Khat / pi)
  LminusR <- Lhat - Kest_obj$r

  data.frame(
    r = Kest_obj$r,
    K = Khat,
    L_minus_r = LminusR
  )
}

# Theoretical curves
k_theory_df <- bind_rows(
  data.frame(
    model = "LGCP",
    r = r_vals,
    K = k_from_pcf(r_vals, function(x) pcf_lgcp_theoretical(x, phi_target, s_target_LGCP))
  ),
  data.frame(
    model = "CSCP",
    r = r_vals,
    K = k_from_pcf(r_vals, function(x) pcf_cscp_theoretical(x, phi_target, s_target))
  )
) %>%
  mutate(
    K_pois = pi * r^2,
    L_minus_r = sqrt(K / pi) - r
  )

# Simulations
results <- vector("list", length = 2 * n_rep)
counter <- 1

for (i in seq_len(n_rep)) {
  cat("Simulation", i, "of", n_rep, "\n")

  X_lgcp <- sim_lgcp_matched(W, lambda_target, phi_target, s_target_LGCP)
  X_cscp <- sim_cscp_matched(W, lambda_target, phi_target, s_target)

  k_lgcp <- estimate_k_curve(X_lgcp, r_vals, correction = "isotropic") %>%
    mutate(model = "LGCP", rep = i)

  k_cscp <- estimate_k_curve(X_cscp, r_vals, correction = "isotropic") %>%
    mutate(model = "CSCP", rep = i)

  results[[counter]] <- k_lgcp
  counter <- counter + 1

  results[[counter]] <- k_cscp
  counter <- counter + 1
}

results_df <- bind_rows(results)

# Aggregate empirical curves
summary_df <- results_df %>%
  group_by(model, r) %>%
  summarise(
    mean_K = mean(K, na.rm = TRUE),
    sd_K   = sd(K, na.rm = TRUE),
    mean_L_minus_r = mean(L_minus_r, na.rm = TRUE),
    sd_L_minus_r   = sd(L_minus_r, na.rm = TRUE),
    .groups = "drop"
  )

print(head(summary_df, 10))

# Plot theoretical K curves
p_theory_K <- ggplot(k_theory_df, aes(x = r, y = K, colour = model)) +
  geom_line(linewidth = 1) +
  geom_line(aes(y = K_pois), linetype = "dashed", colour = "black") +
  labs(
    title = "Theoretical K-functions",
    subtitle = "Dashed line = Poisson benchmark",
    x = "r",
    y = "K(r)"
  ) +
  theme_minimal()

p_theory_L <- ggplot(k_theory_df, aes(x = r, y = L_minus_r, colour = model)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  labs(
    title = "Theoretical L(r) - r curves",
    subtitle = "Zero line = Poisson benchmark",
    x = "r",
    y = "L(r) - r"
  ) +
  theme_minimal()

print(p_theory_K)
print(p_theory_L)

# Plot empirical mean curves
p_emp_K <- ggplot(summary_df, aes(x = r, y = mean_K, colour = model)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Mean empirical K-functions",
    subtitle = paste0("n_rep = ", n_rep),
    x = "r",
    y = "Mean K(r)"
  ) +
  theme_minimal()

p_emp_L <- ggplot(summary_df, aes(x = r, y = mean_L_minus_r, colour = model)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Mean empirical L(r) - r curves",
    subtitle = paste0("n_rep = ", n_rep),
    x = "r",
    y = "Mean L(r) - r"
  ) +
  theme_minimal()

print(p_emp_K)
print(p_emp_L)

# Overlay theoretical and empirical
plot_df_K <- summary_df %>%
  select(model, r, empirical = mean_K) %>%
  left_join(
    k_theory_df %>% select(model, r, theoretical = K),
    by = c("model", "r")
  ) %>%
  pivot_longer(cols = c(empirical, theoretical),
               names_to = "type", values_to = "K")

plot_df_L <- summary_df %>%
  select(model, r, empirical = mean_L_minus_r) %>%
  left_join(
    k_theory_df %>% select(model, r, theoretical = L_minus_r),
    by = c("model", "r")
  ) %>%
  pivot_longer(cols = c(empirical, theoretical),
               names_to = "type", values_to = "L_minus_r")

p_overlay_K <- ggplot(plot_df_K, aes(x = r, y = K, colour = model, linetype = type)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Empirical vs theoretical K-functions",
    x = "r",
    y = "K(r)"
  ) +
  theme_minimal()

p_overlay_L <- ggplot(plot_df_L, aes(x = r, y = L_minus_r, colour = model, linetype = type)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Empirical vs theoretical L(r) - r curves",
    x = "r",
    y = "L(r) - r"
  ) +
  theme_minimal()

print(p_overlay_K)
print(p_overlay_L)

# Simple curve-distance summaries
curve_diff_df <- summary_df %>%
  select(model, r, mean_K, mean_L_minus_r) %>%
  pivot_wider(names_from = model, values_from = c(mean_K, mean_L_minus_r)) %>%
  mutate(
    abs_diff_K = abs(mean_K_LGCP - mean_K_CSCP),
    abs_diff_L = abs(mean_L_minus_r_LGCP - mean_L_minus_r_CSCP)
  )

# Approximate integrated absolute difference using grid spacing
dr <- diff(r_vals)[1]

distance_summary <- data.frame(
  integrated_abs_diff_K = sum(curve_diff_df$abs_diff_K) * dr,
  integrated_abs_diff_L = sum(curve_diff_df$abs_diff_L) * dr,
  max_abs_diff_K = max(curve_diff_df$abs_diff_K),
  max_abs_diff_L = max(curve_diff_df$abs_diff_L)
)

print(distance_summary)

# Optional: inspect one realisation
X_lgcp_demo <- sim_lgcp_matched(W, lambda_target, phi_target, s_target_LGCP)
X_cscp_demo <- sim_cscp_matched(W, lambda_target, phi_target, s_target)

par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
plot(X_cscp_demo, main = "CSCP example")
plot(X_lgcp_demo, main = "LGCP example")
