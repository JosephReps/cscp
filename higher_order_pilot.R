# Pilot: extreme-region geometry from smoothed intensity
rm(list = ls())

library(spatstat)
library(dplyr)
library(ggplot2)
library(tidyr)
devtools::load_all()

set.seed(1)

# Inputs
W <- spatstat.geom::owin(c(0, 1), c(0, 1))

lambda_target <- 300
phi_target    <- 1.5
s_target      <- 0.05
s_target_LGCP <- s_target * 0.68

n_rep <- 200

# KDE bandwidth: start around the CSCP scale
sigma_kde <- 0.04

# Threshold quantiles for excursion sets
thr_probs <- c(0.95, 0.99)

simulate_one_lgcp <- function() {
  sim_lgcp(
    W = W,
    lambda = lambda_target,
    phi = phi_target,
    s = s_target_LGCP
  )$pp
}

simulate_one_cscp <- function() {
  sim_cscp(
    W = W,
    lambda = lambda_target,
    phi = phi_target,
    s = s_target,
    df = 1
  )$pp
}

# Helpers

# simple Gini coefficient
gini_coef <- function(x) {
  x <- x[is.finite(x)]
  x <- x[x >= 0]
  if (length(x) == 0 || mean(x) == 0) return(NA_real_)
  x <- sort(x)
  n <- length(x)
  (2 * sum(seq_len(n) * x)) / (n * sum(x)) - (n + 1) / n
}

# rough skewness proxy
skew3 <- function(x) {
  x <- x[is.finite(x)]
  s <- sd(x)
  if (length(x) < 3 || s == 0) return(NA_real_)
  mean((x - mean(x))^3) / s^3
}

# connected components in a binary matrix using 8-neighbourhood
label_components <- function(mat) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  lab <- matrix(0L, nr, nc)
  current <- 0L

  nbrs <- function(i, j) {
    ii <- pmax(1, i - 1):pmin(nr, i + 1)
    jj <- pmax(1, j - 1):pmin(nc, j + 1)
    expand.grid(ii = ii, jj = jj)
  }

  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      if (!mat[i, j] || lab[i, j] != 0L) next
      current <- current + 1L
      q_i <- i
      q_j <- j
      head <- 1L

      lab[i, j] <- current

      while (head <= length(q_i)) {
        ci <- q_i[head]
        cj <- q_j[head]
        head <- head + 1L

        nn <- nbrs(ci, cj)
        for (k in seq_len(nrow(nn))) {
          ni <- nn$ii[k]
          nj <- nn$jj[k]
          if (mat[ni, nj] && lab[ni, nj] == 0L) {
            lab[ni, nj] <- current
            q_i <- c(q_i, ni)
            q_j <- c(q_j, nj)
          }
        }
      }
    }
  }

  lab
}

# extract summary features from smoothed intensity image
extract_surface_features <- function(X, sigma_kde, thr_probs = c(0.95, 0.99)) {
  lam_hat <- density.ppp(
    X,
    sigma = sigma_kde,
    at = "pixels",
    edge = TRUE
  )

  vals <- as.vector(lam_hat$v)
  vals <- vals[is.finite(vals)]

  out <- tibble(
    mean_hat = mean(vals),
    sd_hat   = sd(vals),
    cv_hat   = sd(vals) / mean(vals),
    gini_hat = gini_coef(vals),
    skew_hat = skew3(vals),
    q95_hat  = unname(quantile(vals, 0.95)),
    q99_hat  = unname(quantile(vals, 0.99)),
    max_hat  = max(vals),
    top1_mass_share = sum(vals[vals >= quantile(vals, 0.99)]) / sum(vals)
  )

  # excursion set features
  for (p in thr_probs) {
    thr <- unname(quantile(vals, p))
    hot <- lam_hat$v >= thr
    hot[!is.finite(hot)] <- FALSE

    labs <- label_components(hot)
    comp_ids <- sort(unique(as.vector(labs)))
    comp_ids <- comp_ids[comp_ids > 0]

    if (length(comp_ids) == 0) {
      n_comp <- 0
      largest_area <- 0
      mean_area <- 0
    } else {
      areas <- sapply(comp_ids, function(id) sum(labs == id))
      n_comp <- length(areas)
      largest_area <- max(areas)
      mean_area <- mean(areas)
    }

    out[[paste0("n_comp_q", 100 * p)]] <- n_comp
    out[[paste0("largest_comp_q", 100 * p)]] <- largest_area
    out[[paste0("mean_comp_q", 100 * p)]] <- mean_area
  }

  out
}

# Run pilot
results <- vector("list", 2 * n_rep)
idx <- 1L

for (i in seq_len(n_rep)) {
  if (i %% 25 == 0) message("rep ", i, " / ", n_rep)

  X_cscp <- simulate_one_cscp()
  X_lgcp <- simulate_one_lgcp()

  results[[idx]] <- extract_surface_features(X_cscp, sigma_kde, thr_probs) %>%
    mutate(model = "CSCP", rep = i)
  idx <- idx + 1L

  results[[idx]] <- extract_surface_features(X_lgcp, sigma_kde, thr_probs) %>%
    mutate(model = "LGCP", rep = i)
  idx <- idx + 1L
}

feat_df <- bind_rows(results)

# Summary table

summary_df <- feat_df %>%
  group_by(model) %>%
  summarise(across(
    .cols = -rep,
    .fns = mean,
    .names = "mean_{.col}"
  ), .groups = "drop")

print(summary_df)

# Pair by rep

pair_df <- feat_df %>%
  pivot_wider(
    id_cols = rep,
    names_from = model,
    values_from = -c(model, rep)
  )

compare_prob <- function(a, b) mean(a > b, na.rm = TRUE)

classifier_summary <- tibble(
  stat = c(
    "cv_hat", "gini_hat", "skew_hat", "q95_hat", "q99_hat", "max_hat",
    "top1_mass_share", "n_comp_q95", "largest_comp_q95",
    "n_comp_q99", "largest_comp_q99"
  ),
  lgcp_gt_cscp = c(
    compare_prob(pair_df$cv_hat_LGCP, pair_df$cv_hat_CSCP),
    compare_prob(pair_df$gini_hat_LGCP, pair_df$gini_hat_CSCP),
    compare_prob(pair_df$skew_hat_LGCP, pair_df$skew_hat_CSCP),
    compare_prob(pair_df$q95_hat_LGCP, pair_df$q95_hat_CSCP),
    compare_prob(pair_df$q99_hat_LGCP, pair_df$q99_hat_CSCP),
    compare_prob(pair_df$max_hat_LGCP, pair_df$max_hat_CSCP),
    compare_prob(pair_df$top1_mass_share_LGCP, pair_df$top1_mass_share_CSCP),
    compare_prob(pair_df$n_comp_q95_LGCP, pair_df$n_comp_q95_CSCP),
    compare_prob(pair_df$largest_comp_q95_LGCP, pair_df$largest_comp_q95_CSCP),
    compare_prob(pair_df$n_comp_q99_LGCP, pair_df$n_comp_q99_CSCP),
    compare_prob(pair_df$largest_comp_q99_LGCP, pair_df$largest_comp_q99_CSCP)
  )
)

print(classifier_summary)

# Boxplots


plot_stats <- c(
  "cv_hat", "gini_hat", "skew_hat", "q95_hat", "q99_hat", "max_hat",
  "top1_mass_share", "n_comp_q95", "largest_comp_q95",
  "n_comp_q99", "largest_comp_q99"
)

plot_df <- feat_df %>%
  select(model, rep, all_of(plot_stats)) %>%
  pivot_longer(
    cols = all_of(plot_stats),
    names_to = "stat",
    values_to = "value"
  )

p_box <- ggplot(plot_df, aes(x = model, y = value, fill = model)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~ stat, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    title = "Smoothed intensity surface features",
    subtitle = paste0("sigma_kde = ", sigma_kde),
    x = NULL,
    y = NULL
  )

print(p_box)

# one example pair

X_cscp_demo <- simulate_one_cscp()
X_lgcp_demo <- simulate_one_lgcp()

lam_cscp <- density.ppp(X_cscp_demo, sigma = sigma_kde, at = "pixels", edge = TRUE)
lam_lgcp <- density.ppp(X_lgcp_demo, sigma = sigma_kde, at = "pixels", edge = TRUE)

par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
plot(X_cscp_demo, main = "CSCP pattern")
plot(X_lgcp_demo, main = "LGCP pattern")
plot(lam_cscp, main = "CSCP KDE intensity")
plot(lam_lgcp, main = "LGCP KDE intensity")
