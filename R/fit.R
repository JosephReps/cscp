## Pending documentation
# See thesis notes for more information

# Pair correlation estimation helper
estimate_pcf_for_mc <- function(X,
                                bw = 0.2,
                                correction = "isotropic",
                                divisor = "a",
                                zerocor = "JonesFoster",
                                r_min = 0.02,
                                r_max = NULL) {
  stopifnot(inherits(X, "ppp"))

  obs <- spatstat.explore::pcf(
    X,
    bw = bw,
    correction = correction,
    divisor = divisor,
    zerocor = zerocor
  )

  y_name <- switch(
    correction,
    isotropic = "iso",
    translate = "trans",
    Ripley = "iso",
    "iso"
  )

  if (!y_name %in% names(obs)) {
    y_name <- intersect(c("iso", "trans", "un"), names(obs))[1]
  }

  r <- obs$r
  g <- obs[[y_name]]

  keep <- is.finite(r) & is.finite(g) & r > r_min
  if (!is.null(r_max)) keep <- keep & r <= r_max

  list(
    fv = obs,
    r = r[keep],
    g = g[keep],
    correction_used = y_name,
    bw = bw,
    divisor = divisor,
    zerocor = zerocor,
    r_min = r_min,
    r_max = r_max
  )
}

# Minimum contrast objective on the PCF scale
mc_pcf_objective <- function(par,
                             r,
                             g_obs,
                             model = c("lgcp", "cscp"),
                             lambda,
                             df = 1,
                             k = 0.25,
                             p = 2) {
  model <- match.arg(model)

  if (model == "lgcp") {
    phi <- exp(par[1])              # phi in (0, Inf)
    s   <- exp(par[2])              # s > 0
    g_th <- pcf_lgcp_theoretical(
      r,
      lambda = lambda,
      phi = phi,
      s = s
    )
  } else {
    eta <- par[1]                   # eta in R
    phi <- 2 * stats::plogis(eta)   # phi in (0, 2)
    s   <- exp(par[2])              # s > 0
    g_th <- pcf_cscp_theoretical(
      r,
      lambda = lambda,
      phi = phi,
      s = s,
      df = df
    )
  }

  sum(abs(g_obs^k - g_th^k)^p)
}

# MC fitting wrappers

fit_lgcp_mc <- function(X,
                        start = list(phi = 0.5, s = 0.05),
                        bw = 0.2,
                        correction = "isotropic",
                        divisor = "a",
                        zerocor = "JonesFoster",
                        r_min = 0.02,
                        r_max = NULL,
                        k = 0.25,
                        p = 2) {
  stopifnot(inherits(X, "ppp"))
  stopifnot(is.list(start), !is.null(start$phi), !is.null(start$s))
  stopifnot(start$phi > 0, start$s > 0)

  lambda_hat <- spatstat.geom::intensity(X)

  pcf_obs <- estimate_pcf_for_mc(
    X,
    bw = bw,
    correction = correction,
    divisor = divisor,
    zerocor = zerocor,
    r_min = r_min,
    r_max = r_max
  )

  opt <- stats::optim(
    par = c(log(start$phi), log(start$s)),
    fn = mc_pcf_objective,
    r = pcf_obs$r,
    g_obs = pcf_obs$g,
    model = "lgcp",
    lambda = lambda_hat,
    k = k,
    p = p,
    method = "L-BFGS-B",
    lower = c(log(1e-8), log(1e-8)),
    upper = c(log(1e3), log(10))
  )

  pars <- c(
    lambda = lambda_hat,
    phi = exp(opt$par[1]),
    s = exp(opt$par[2])
  )

  g_fit <- pcf_lgcp_theoretical(
    pcf_obs$r,
    lambda = pars["lambda"],
    phi = pars["phi"],
    s = pars["s"]
  )

  structure(
    list(
      model = "LGCP",
      par = pars,
      objective = opt$value,
      convergence = opt$convergence,
      optim = opt,
      r = pcf_obs$r,
      g_obs = pcf_obs$g,
      g_fit = g_fit,
      pcf_obs = pcf_obs
    ),
    class = "mc_pcf_fit"
  )
}

fit_cscp_mc <- function(X,
                        start = list(phi = 0.5, s = 0.05),
                        df = 1,
                        bw = 0.2,
                        correction = "isotropic",
                        divisor = "a",
                        zerocor = "JonesFoster",
                        r_min = 0.02,
                        r_max = NULL,
                        k = 0.25,
                        p = 2) {
  stopifnot(inherits(X, "ppp"))
  stopifnot(is.list(start), !is.null(start$phi), !is.null(start$s))
  stopifnot(start$phi > 0, start$phi < 2, start$s > 0)

  lambda_hat <- spatstat.geom::intensity(X)

  pcf_obs <- estimate_pcf_for_mc(
    X,
    bw = bw,
    correction = correction,
    divisor = divisor,
    zerocor = zerocor,
    r_min = r_min,
    r_max = r_max
  )

  # Inverse of phi = 2*plogis(eta)
  eta_start <- stats::qlogis(start$phi / 2)

  opt <- stats::optim(
    par = c(eta_start, log(start$s)),
    fn = mc_pcf_objective,
    r = pcf_obs$r,
    g_obs = pcf_obs$g,
    model = "cscp",
    lambda = lambda_hat,
    df = df,
    k = k,
    p = p,
    method = "L-BFGS-B",
    lower = c(-20, log(1e-8)),
    upper = c(20, log(10))
  )

  pars <- c(
    lambda = lambda_hat,
    phi = 2 * stats::plogis(opt$par[1]),
    s = exp(opt$par[2])
  )

  g_fit <- pcf_cscp_theoretical(
    pcf_obs$r,
    lambda = pars["lambda"],
    phi = pars["phi"],
    s = pars["s"],
    df = df
  )

  structure(
    list(
      model = "CSCP",
      df = df,
      par = pars,
      objective = opt$value,
      convergence = opt$convergence,
      optim = opt,
      eta_hat = opt$par[1],
      r = pcf_obs$r,
      g_obs = pcf_obs$g,
      g_fit = g_fit,
      pcf_obs = pcf_obs
    ),
    class = "mc_pcf_fit"
  )
}
