# Pilot: synthetic-likelihood discrimination

rm(list = ls())
library(spatstat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mvtnorm)
library(MASS)

set.seed(1)

# Inputs
W <- spatstat.geom::owin(c(0, 1), c(0, 1))

lambda_target <- 100
phi_target    <- 1.5
s_target      <- 0.05
s_target_LGCP <- s_target * 0.68

# Number of test datasets to classify
n_test <- 10

# Number of simulations used to estimate synthetic likelihood
n_sim_per_candidate <- 30

# Grid resolution for count field
nx_grid <- 12
ny_grid <- 12

# Candidate parameter grid
phi_grid <- c(0.75, 1.0, 1.5, 2.0)
s_grid_cscp <- c(0.03, 0.05, 0.08)
s_grid_lgcp <- s_grid_cscp * 0.68

# Small ridge for covariance stability
ridge_eps <- 1e-6

# Wrappers
# Replace these two if your sim function names differ
simulate_one_lgcp <- function(lambda, phi, s_lgcp) {
  sim_lgcp(
    W = W,
    lambda = lambda,
    phi = phi,
    s = s_lgcp
  )$pp
}

simulate_one_cscp <- function(lambda, phi, s_cscp) {
  sim_cscp(
    W = W,
    lambda = lambda,
    phi = phi,
    s = s_cscp,
    df = 1
  )$pp
}

# Helpers

# Convert point pattern to counts on fixed grid
pattern_to_count_matrix <- function(X, nx = nx_grid, ny = ny_grid) {
  qc <- quadratcount(X, nx = nx, ny = ny)
  mat <- matrix(as.numeric(qc), nrow = ny, ncol = nx, byrow = TRUE)
  mat
}

# Aggregate count matrix into 2x2 blocks
aggregate_2x2 <- function(mat) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  nr2 <- floor(nr / 2)
  nc2 <- floor(nc / 2)

  out <- matrix(0, nrow = nr2, ncol = nc2)

  for (i in seq_len(nr2)) {
    for (j in seq_len(nc2)) {
      rows <- (2 * i - 1):(2 * i)
      cols <- (2 * j - 1):(2 * j)
      out[i, j] <- sum(mat[rows, cols])
    }
  }
  out
}

# Simple skewness
skew3 <- function(x) {
  s <- sd(x)
  if (length(x) < 3 || s == 0) return(0)
  mean((x - mean(x))^3) / s^3
}

# Horizontal and vertical nearest-neighbour correlation
adjacent_corr <- function(mat) {
  h_x <- as.vector(mat[, -ncol(mat)])
  h_y <- as.vector(mat[, -1])

  v_x <- as.vector(mat[-nrow(mat), ])
  v_y <- as.vector(mat[-1, ])

  h_cor <- suppressWarnings(cor(h_x, h_y))
  v_cor <- suppressWarnings(cor(v_x, v_y))

  h_cor <- ifelse(is.finite(h_cor), h_cor, 0)
  v_cor <- ifelse(is.finite(v_cor), v_cor, 0)

  c(h_cor = h_cor, v_cor = v_cor)
}

# Summary vector from one pattern
# This is the key object for synthetic likelihood
pattern_summary <- function(X) {
  mat1 <- pattern_to_count_matrix(X, nx_grid, ny_grid)
  vec1 <- as.vector(mat1)

  mat2 <- aggregate_2x2(mat1)
  vec2 <- as.vector(mat2)

  adj <- adjacent_corr(mat1)

  c(
    mean1      = mean(vec1),
    var1       = var(vec1),
    skew1      = skew3(vec1),
    p0_1       = mean(vec1 == 0),
    pge3_1     = mean(vec1 >= 3),
    pge5_1     = mean(vec1 >= 5),
    q95_1      = unname(quantile(vec1, 0.95)),
    max1       = max(vec1),

    mean2      = mean(vec2),
    var2       = var(vec2),
    skew2      = skew3(vec2),
    p0_2       = mean(vec2 == 0),
    pge8_2     = mean(vec2 >= 8),
    pge12_2    = mean(vec2 >= 12),
    q95_2      = unname(quantile(vec2, 0.95)),
    max2       = max(vec2),

    h_cor      = adj["h_cor"],
    v_cor      = adj["v_cor"]
  )
}

# Simulate summary vectors under candidate model
simulate_summary_matrix <- function(model, lambda, phi, s, n_sim) {
  out <- vector("list", n_sim)

  for (i in seq_len(n_sim)) {
    X <- if (model == "LGCP") {
      simulate_one_lgcp(lambda = lambda, phi = phi, s_lgcp = s)
    } else {
      simulate_one_cscp(lambda = lambda, phi = phi, s_cscp = s)
    }

    out[[i]] <- pattern_summary(X)
  }

  out <- do.call(rbind, out)
  out
}

# Synthetic log-likelihood:
# Assume summaries ~ approximately MVN under candidate
synthetic_loglik <- function(s_obs, S_sim, ridge = ridge_eps) {
  mu <- colMeans(S_sim)
  Sigma <- stats::cov(S_sim)

  # ridge stabilize
  Sigma <- Sigma + diag(ridge, ncol(Sigma))

  ll <- tryCatch(
    dmvnorm(s_obs, mean = mu, sigma = Sigma, log = TRUE),
    error = function(e) NA_real_
  )

  ll
}

# Fit candidate family by grid search
fit_family_synthlik <- function(X, family = c("LGCP", "CSCP")) {
  family <- match.arg(family)

  # Estimate lambda from observed count
  lambda_hat <- npoints(X) / area.owin(W)

  s_obs <- pattern_summary(X)

  if (family == "LGCP") {
    grid <- expand.grid(phi = phi_grid, s = s_grid_lgcp)
  } else {
    grid <- expand.grid(phi = phi_grid, s = s_grid_cscp)
  }

  scores <- numeric(nrow(grid))

  for (k in seq_len(nrow(grid))) {
    phi_k <- grid$phi[k]
    s_k   <- grid$s[k]

    S_sim <- simulate_summary_matrix(
      model = family,
      lambda = lambda_hat,
      phi = phi_k,
      s = s_k,
      n_sim = n_sim_per_candidate
    )

    scores[k] <- synthetic_loglik(s_obs, S_sim)
  }

  best_idx <- which.max(scores)

  list(
    family = family,
    lambda_hat = lambda_hat,
    s_obs = s_obs,
    grid = grid,
    scores = scores,
    best_phi = grid$phi[best_idx],
    best_s = grid$s[best_idx],
    best_score = scores[best_idx]
  )
}

# 1) Single-dataset pilot
cat("---- SINGLE DATASET PILOT ----\n")

X_true <- simulate_one_lgcp(
  lambda = lambda_target,
  phi = phi_target,
  s_lgcp = s_target_LGCP
)

fit_lgcp <- fit_family_synthlik(X_true, family = "LGCP")
fit_cscp <- fit_family_synthlik(X_true, family = "CSCP")

single_result <- tibble(
  true_model = "LGCP",
  fit_family = c("LGCP", "CSCP"),
  best_phi   = c(fit_lgcp$best_phi, fit_cscp$best_phi),
  best_s     = c(fit_lgcp$best_s, fit_cscp$best_s),
  best_score = c(fit_lgcp$best_score, fit_cscp$best_score)
)

print(single_result)

# 2) Repeated classification study
cat("---- REPEATED CLASSIFICATION STUDY ----\n")

results <- vector("list", 2 * n_test)
idx <- 1L

for (i in seq_len(n_test)) {
  cat("Test dataset", i, "of", n_test, "\n")

  # -------- true LGCP dataset --------
  X_lgcp <- simulate_one_lgcp(
    lambda = lambda_target,
    phi = phi_target,
    s_lgcp = s_target_LGCP
  )

  fit_l_on_l <- fit_family_synthlik(X_lgcp, family = "LGCP")
  fit_c_on_l <- fit_family_synthlik(X_lgcp, family = "CSCP")

  results[[idx]] <- tibble(
    rep = i,
    true_model = "LGCP",
    score_LGCP = fit_l_on_l$best_score,
    score_CSCP = fit_c_on_l$best_score,
    chosen_model = ifelse(fit_l_on_l$best_score > fit_c_on_l$best_score, "LGCP", "CSCP"),
    best_phi_LGCP = fit_l_on_l$best_phi,
    best_s_LGCP = fit_l_on_l$best_s,
    best_phi_CSCP = fit_c_on_l$best_phi,
    best_s_CSCP = fit_c_on_l$best_s
  )
  idx <- idx + 1L

  # -------- true CSCP dataset --------
  X_cscp <- simulate_one_cscp(
    lambda = lambda_target,
    phi = phi_target,
    s_cscp = s_target
  )

  fit_l_on_c <- fit_family_synthlik(X_cscp, family = "LGCP")
  fit_c_on_c <- fit_family_synthlik(X_cscp, family = "CSCP")

  results[[idx]] <- tibble(
    rep = i,
    true_model = "CSCP",
    score_LGCP = fit_l_on_c$best_score,
    score_CSCP = fit_c_on_c$best_score,
    chosen_model = ifelse(fit_l_on_c$best_score > fit_c_on_c$best_score, "LGCP", "CSCP"),
    best_phi_LGCP = fit_l_on_c$best_phi,
    best_s_LGCP = fit_l_on_c$best_s,
    best_phi_CSCP = fit_c_on_c$best_phi,
    best_s_CSCP = fit_c_on_c$best_s
  )
  idx <- idx + 1L
}

results_df <- bind_rows(results) %>%
  mutate(correct = chosen_model == true_model,
         score_diff = score_LGCP - score_CSCP)

print(results_df)

# 3) Summary tables
summary_table <- results_df %>%
  group_by(true_model) %>%
  summarise(
    n = n(),
    accuracy = mean(correct, na.rm = TRUE),
    mean_score_diff = mean(score_diff, na.rm = TRUE),
    median_score_diff = median(score_diff, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_table)

overall_accuracy <- mean(results_df$correct, na.rm = TRUE)
cat("\nOverall classification accuracy:", round(overall_accuracy, 3), "\n")

# 4) Plots
p_scores <- ggplot(results_df, aes(x = true_model, y = score_diff, fill = true_model)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Score difference: LGCP synthetic likelihood - CSCP synthetic likelihood",
    subtitle = "Positive values favour LGCP; negative values favour CSCP",
    x = "True generating model",
    y = "Score difference"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p_scores)

p_choice <- ggplot(results_df, aes(x = true_model, fill = chosen_model)) +
  geom_bar(position = "fill") +
  labs(
    title = "Chosen model by true generating model",
    x = "True generating model",
    y = "Proportion"
  ) +
  theme_minimal()

print(p_choice)
