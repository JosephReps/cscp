devtools::load_all()
library(spatstat)

simulate_one_model <- function(model, W, lambda, phi, s_cscp, df = 1) {
  if (model == "CSCP") {
    out <- sim_cscp(
      W = W,
      lambda = lambda,
      phi = phi,
      s = s_cscp,
      df = df
    )
  } else if (model == "LGCP") {
    s_lgcp <- match_lgcp_scale(phi = phi, s_cscp = s_cscp)
    out <- sim_lgcp(
      W = W,
      lambda = lambda,
      phi = phi,
      s = s_lgcp
    )
  } else {
    stop("model must be 'LGCP' or 'CSCP'")
  }

  out
}

make_shuffled_labels <- function(n_each = 8) {
  sample(rep(c("LGCP", "CSCP"), each = n_each))
}

generate_regime_draws <- function(W, lambda, phi, s_cscp, df = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  labels_surface <- make_shuffled_labels(8)
  labels_pattern <- make_shuffled_labels(8)

  surface_sims <- lapply(labels_surface, function(m) {
    simulate_one_model(
      model = m, W = W, lambda = lambda, phi = phi,
      s_cscp = s_cscp, df = df
    )
  })

  pattern_sims <- lapply(labels_pattern, function(m) {
    simulate_one_model(
      model = m, W = W, lambda = lambda, phi = phi,
      s_cscp = s_cscp, df = df
    )
  })

  list(
    settings = list(lambda = lambda, phi = phi, s_cscp = s_cscp, df = df),
    surface_labels = labels_surface,
    pattern_labels = labels_pattern,
    surface_sims = surface_sims,
    pattern_sims = pattern_sims
  )
}

get_shared_zlim <- function(sim_list) {
  vals <- unlist(lapply(sim_list, function(x) as.vector(x$intensity$v)))
  range(vals, finite = TRUE)
}

plot_surface_grid <- function(sim_list,
                              col_fun = NULL,
                              zlim = NULL,
                              mar = c(0.1, 0.1, 0.1, 0.1),
                              oma = c(0, 0, 0, 0)) {

  if (is.null(col_fun)) {
    # close to spatstat-style heat colours
    col_fun <- function(n) {
      rev(heat.colors(n))
    }
  }

  if (is.null(zlim)) {
    zlim <- get_shared_zlim(sim_list)
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow = c(4, 4), mar = mar, oma = oma)

  for (i in seq_along(sim_list)) {
    plot(
      sim_list[[i]]$intensity,
      ribbon = FALSE,
      zlim = zlim,
      axes = FALSE,
      box = FALSE,
      main = "",
      ann = FALSE
    )
  }

  invisible(zlim)
}

plot_pattern_grid <- function(sim_list,
                              pch = 16,
                              cex = 0.35,
                              mar = c(0.1, 0.1, 0.1, 0.1),
                              oma = c(0, 0, 0, 0)) {

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow = c(4, 4), mar = mar, oma = oma)

  for (i in seq_along(sim_list)) {
    plot(
      sim_list[[i]]$pp,
      main = "",
      axes = FALSE,
      ann = FALSE,
      chars = pch,
      cols = "black",
      cex = cex
    )
    box()
  }
}

print_answer_key <- function(regime_obj, regime_name = "Regime") {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat(regime_name, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat("Parameters:\n")
  cat(sprintf(
    "  lambda = %s, phi = %s, s_cscp = %s, df = %s\n\n",
    regime_obj$settings$lambda,
    regime_obj$settings$phi,
    regime_obj$settings$s_cscp,
    regime_obj$settings$df
  ))

  surf_mat <- matrix(regime_obj$surface_labels, nrow = 4, byrow = TRUE)
  patt_mat <- matrix(regime_obj$pattern_labels, nrow = 4, byrow = TRUE)

  cat("Surface grid answer key (row-wise 4x4):\n")
  print(surf_mat, quote = FALSE)

  cat("\nPattern grid answer key (row-wise 4x4):\n")
  print(patt_mat, quote = FALSE)

  invisible(NULL)
}

make_guessing_game_grids <- function(
    regimes = list(
      list(lambda = 100,  phi = 0.5, s_cscp = 0.05),
      list(lambda = 500,  phi = 1.0, s_cscp = 0.05),
      list(lambda = 2000, phi = 2.0, s_cscp = 0.05)
    ),
    W = spatstat.geom::owin(c(0, 1), c(0, 1)),
    df = 1,
    seed = 123,
    surface_cex_comment = FALSE
) {
  out <- vector("list", length(regimes))

  for (r in seq_along(regimes)) {
    reg <- regimes[[r]]

    out[[r]] <- generate_regime_draws(
      W = W,
      lambda = reg$lambda,
      phi = reg$phi,
      s_cscp = reg$s_cscp,
      df = df,
      seed = seed + 1000 * r
    )
  }

  names(out) <- paste0("Regime_", seq_along(regimes))
  out
}

save_regime_plots <- function(regime_obj,
                              prefix = "regime1",
                              width = 8,
                              height = 8,
                              dpi = 300) {

  zlim <- get_shared_zlim(regime_obj$surface_sims)

  png(
    filename = paste0(prefix, "_surfaces.png"),
    width = width,
    height = height,
    units = "in",
    res = dpi
  )
  plot_surface_grid(
    regime_obj$surface_sims,
    zlim = zlim
  )
  dev.off()

  png(
    filename = paste0(prefix, "_patterns.png"),
    width = width,
    height = height,
    units = "in",
    res = dpi
  )
  plot_pattern_grid(
    regime_obj$pattern_sims
  )
  dev.off()
}

regimes <- list(
  list(lambda = 500,  phi = 0.5, s_cscp = 0.05),
  list(lambda = 500,  phi = 1.0, s_cscp = 0.12),
  list(lambda = 500, phi = 2.0, s_cscp = 0.22)
)

game <- make_guessing_game_grids(
  regimes = regimes,
  W = spatstat.geom::owin(c(0, 1), c(0, 1)),
  df = 1,
  seed = 2026
)

# Plot them to screen
# For each regime:
#   1) surfaces
#   2) patterns
# Shared colour scale is enforced within each regime's 16 surfaces.

for (r in seq_along(game)) {
  regime_name <- names(game)[r]
  cat("\nShowing", regime_name, "...\n")

  zlim_regime <- get_shared_zlim(game[[r]]$surface_sims)

  # surface grid
  # plot_surface_grid(
  #   game[[r]]$surface_sims,
  #   zlim = zlim_regime
  # )

  # point pattern grid
  plot_pattern_grid(
    game[[r]]$pattern_sims
  )

  # print answer key
  print_answer_key(
    game[[r]],
    regime_name = regime_name
  )
}

# save all to PNGs
# for (r in seq_along(game)) {
#   save_regime_plots(
#     regime_obj = game[[r]],
#     prefix = paste0("guess_game_", names(game)[r])
#   )
# }
