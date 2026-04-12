# Pilot: distribution of quadrat counts
rm(list = ls())

library(spatstat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
devtools::load_all()

set.seed(1)

# Window + matched parameters
W <- spatstat.geom::owin(c(0, 1), c(0, 1))

lambda_target <- 1000
phi_target    <- 1.5
s_target      <- 0.05
s_target_LGCP <- s_target * 0.68   # IMPORTANT: matched LGCP scale

n_rep <- 200

# Main quadrat size for the first pilot
# Chosen to be large enough to see aggregation, but not so large that
# everything gets averaged away.
q_main <- 0.10

# Small sensitivity grid around the main choice
q_grid <- c(0.06, 0.08, 0.10, 0.12, 0.15)

# Helpers

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

# Quadrat counts for side length q
get_quadrat_counts <- function(X, q) {
  xr <- diff(X$window$xrange)
  yr <- diff(X$window$yrange)

  nx <- max(1, floor(xr / q))
  ny <- max(1, floor(yr / q))

  qc <- quadratcount(X, nx = nx, ny = ny)

  tibble(
    count = as.numeric(qc),
    nx = nx,
    ny = ny,
    q = q
  )
}

# Summaries of the full count distribution
summarise_count_dist <- function(counts_vec) {
  tibble(
    mean_count   = mean(counts_vec),
    var_count    = var(counts_vec),
    fano         = var(counts_vec) / mean(counts_vec),
    skewness_proxy = mean((counts_vec - mean(counts_vec))^3) / (sd(counts_vec)^3 + 1e-8),
    q90          = unname(quantile(counts_vec, 0.90)),
    q95          = unname(quantile(counts_vec, 0.95)),
    q99          = unname(quantile(counts_vec, 0.99)),
    max_count    = max(counts_vec),
    p_ge_5       = mean(counts_vec >= 5),
    p_ge_10      = mean(counts_vec >= 10),
    p_ge_15      = mean(counts_vec >= 15),
    p_ge_20      = mean(counts_vec >= 20)
  )
}

# 1) Main pilot at one fixed quadrat size
message("Running main pilot at q = ", q_main)

all_counts_main <- vector("list", 2 * n_rep)
all_summaries_main <- vector("list", 2 * n_rep)
idx <- 1

for (i in seq_len(n_rep)) {
  if (i %% 25 == 0) message("rep ", i, " / ", n_rep)

  X_cscp <- simulate_one_cscp()
  X_lgcp <- simulate_one_lgcp()

  counts_cscp <- get_quadrat_counts(X_cscp, q_main) %>%
    mutate(model = "CSCP", rep = i)

  counts_lgcp <- get_quadrat_counts(X_lgcp, q_main) %>%
    mutate(model = "LGCP", rep = i)

  all_counts_main[[idx]] <- counts_cscp
  idx <- idx + 1
  all_counts_main[[idx]] <- counts_lgcp
  idx <- idx + 1

  all_summaries_main[[idx - 1]] <- summarise_count_dist(counts_cscp$count) %>%
    mutate(model = "LGCP", rep = i)

  all_summaries_main[[idx - 2]] <- summarise_count_dist(counts_lgcp$count) %>%
    mutate(model = "CSCP", rep = i)
}

# fix accidental ordering from above two lines
summary_main <- bind_rows(all_summaries_main) %>%
  mutate(model = ifelse(row_number() %% 2 == 1, "CSCP", "LGCP")) %>%
  group_by(model, rep) %>%
  slice(1) %>%
  ungroup()

counts_main <- bind_rows(all_counts_main)

# 2) Plot one cherry-picked replicate with visible contrast
rep_extreme <- summary_main %>%
  mutate(tail_score = q99 + max_count + p_ge_15 * 100) %>%
  group_by(model) %>%
  slice_max(order_by = tail_score, n = 1) %>%
  ungroup()

rep_lgcp_show <- rep_extreme$rep[rep_extreme$model == "LGCP"][1]
rep_cscp_show <- rep_extreme$rep[rep_extreme$model == "CSCP"][1]

plot_df_show <- counts_main %>%
  filter((model == "LGCP" & rep == rep_lgcp_show) |
           (model == "CSCP" & rep == rep_cscp_show))

p_hist_show <- ggplot(plot_df_show, aes(x = count)) +
  geom_histogram(binwidth = 1, boundary = -0.5, closed = "right") +
  facet_wrap(~ model, ncol = 2, scales = "free_y") +
  labs(
    title = paste0("Quadrat count distribution at q = ", q_main),
    subtitle = "Cherry-picked replicates highlighting tail behaviour",
    x = "Quadrat count",
    y = "Frequency"
  ) +
  theme_minimal()

print(p_hist_show)

# 3) Aggregate histogram across all quadrats from all reps
counts_main_grouped <- counts_main %>%
  group_by(model, count) %>%
  summarise(freq = n(), .groups = "drop") %>%
  group_by(model) %>%
  mutate(prob = freq / sum(freq)) %>%
  ungroup()

p_hist_all <- ggplot(counts_main_grouped, aes(x = count, y = prob, fill = model)) +
  geom_col(position = "identity", alpha = 0.45) +
  labs(
    title = paste0("Pooled quadrat count distribution at q = ", q_main),
    x = "Quadrat count",
    y = "Empirical probability"
  ) +
  theme_minimal()

print(p_hist_all)

# Tail-only version
tail_cut <- quantile(counts_main$count, 0.90)

p_tail_all <- counts_main %>%
  filter(count >= tail_cut) %>%
  ggplot(aes(x = count, fill = model)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.45) +
  labs(
    title = paste0("Upper tail of pooled quadrat counts (q = ", q_main, ")"),
    x = "Quadrat count",
    y = "Frequency"
  ) +
  theme_minimal()

print(p_tail_all)

# 4) Replicate-level summaries
summary_main_by_model <- summary_main %>%
  group_by(model) %>%
  summarise(
    mean_mean_count = mean(mean_count),
    mean_var_count  = mean(var_count),
    mean_fano       = mean(fano),
    mean_skewness   = mean(skewness_proxy),
    mean_q90        = mean(q90),
    mean_q95        = mean(q95),
    mean_q99        = mean(q99),
    mean_max        = mean(max_count),
    mean_p_ge_10    = mean(p_ge_10),
    mean_p_ge_15    = mean(p_ge_15),
    mean_p_ge_20    = mean(p_ge_20),
    .groups = "drop"
  )

print(summary_main_by_model)

# Boxplots of upper-tail summaries across reps
summary_long <- summary_main %>%
  select(model, rep, q95, q99, max_count, p_ge_10, p_ge_15, p_ge_20) %>%
  pivot_longer(
    cols = -c(model, rep),
    names_to = "stat",
    values_to = "value"
  )

p_box <- ggplot(summary_long, aes(x = model, y = value, fill = model)) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(~ stat, scales = "free_y") +
  labs(
    title = paste0("Replicate-level tail summaries at q = ", q_main),
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p_box)

# 5) Tiny classifier pilot

# This is NOT meant to be final inference.
# It just checks whether simple tail statistics carry signal.

classify_stat <- function(x_lgcp, x_cscp) {
  mean(x_lgcp > x_cscp)
}

compare_table <- summary_main %>%
  select(model, rep, q95, q99, max_count, p_ge_10, p_ge_15, p_ge_20) %>%
  pivot_wider(names_from = model, values_from = c(q95, q99, max_count, p_ge_10, p_ge_15, p_ge_20))

classifier_summary <- tibble(
  stat = c("q95", "q99", "max_count", "p_ge_10", "p_ge_15", "p_ge_20"),
  lgcp_gt_cscp = c(
    mean(compare_table$q95_LGCP      > compare_table$q95_CSCP),
    mean(compare_table$q99_LGCP      > compare_table$q99_CSCP),
    mean(compare_table$max_count_LGCP > compare_table$max_count_CSCP),
    mean(compare_table$p_ge_10_LGCP  > compare_table$p_ge_10_CSCP),
    mean(compare_table$p_ge_15_LGCP  > compare_table$p_ge_15_CSCP),
    mean(compare_table$p_ge_20_LGCP  > compare_table$p_ge_20_CSCP)
  )
)

print(classifier_summary)

# 6) Sensitivity to quadrat size
message("Running quadrat-size sensitivity check")

sens_results <- list()
idx <- 1

for (q in q_grid) {
  message("q = ", q)

  for (i in seq_len(n_rep)) {
    X_cscp <- simulate_one_cscp()
    X_lgcp <- simulate_one_lgcp()

    s_cscp <- summarise_count_dist(get_quadrat_counts(X_cscp, q)$count) %>%
      mutate(model = "CSCP", rep = i, q = q)

    s_lgcp <- summarise_count_dist(get_quadrat_counts(X_lgcp, q)$count) %>%
      mutate(model = "LGCP", rep = i, q = q)

    sens_results[[idx]] <- s_cscp
    idx <- idx + 1
    sens_results[[idx]] <- s_lgcp
    idx <- idx + 1
  }
}

sens_df <- bind_rows(sens_results)

sens_summary <- sens_df %>%
  group_by(model, q) %>%
  summarise(
    q95 = mean(q95),
    q99 = mean(q99),
    max_count = mean(max_count),
    p_ge_10 = mean(p_ge_10),
    p_ge_15 = mean(p_ge_15),
    p_ge_20 = mean(p_ge_20),
    .groups = "drop"
  )

print(sens_summary)

sens_long <- sens_summary %>%
  pivot_longer(
    cols = -c(model, q),
    names_to = "stat",
    values_to = "value"
  )

p_sens <- ggplot(sens_long, aes(x = q, y = value, colour = model)) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_wrap(~ stat, scales = "free_y") +
  labs(
    title = "Sensitivity of tail summaries to quadrat size",
    x = "Quadrat side length",
    y = NULL
  ) +
  theme_minimal()

print(p_sens)

# 7) Optional: one example point pattern pair
X_cscp_demo <- simulate_one_cscp()
X_lgcp_demo <- simulate_one_lgcp()

par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
plot(X_cscp_demo, main = "CSCP example")
plot(X_lgcp_demo, main = "LGCP example")


















