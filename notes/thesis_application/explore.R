library(spatstat)
library(dplyr)
library(purrr)

# aegiss <- readRDS("./notes/thesis_application/aegiss_spatiotemporal.rds")
#
# sel <- seq(1, max(marks(aegiss)), by = 6)
# X_sub <- aegiss[marks(aegiss) %in% sel]
# X_sub <- unmark(X_sub)
#
# lambda_hat <- intensity(X)
#
# fit_lgcp <- kppm(
#   X,
#   clusters = "LGCP",
#   statistic = "pcf"
# )
#
# fit_lgcp

X <- bei
bw = 10 # Not passing right now

test_obs <- pcf(X, correction = "isotropic", divisor = "a", zerocor = "JonesFoster", bw = bw)
r <- test_obs$r

test_fit <- fit_lgcp_mc(X, correction = "isotropic", divisor = "a", zerocor = "JonesFoster",
                        bw = bw)
test_fitted_line <- pcf_lgcp_theoretical(r, lambda = test_fit$par[["lambda"]], phi = test_fit$par[["phi"]],
                                         s = test_fit$par[["s"]])
test_fit2 <- fit_cscp_mc(X, correction = "isotropic", divisor = "a", zerocor = "JonesFoster",
                         bw = bw)
test_fitted_line2 <- pcf_lgcp_theoretical(r, lambda = test_fit$par[["lambda"]], phi = test_fit$par[["phi"]],
                                         s = test_fit$par[["s"]])


# Plotting
plot(r, test_obs$iso, type = "l", lwd = 2, lty = "dashed",
     main = "CSCP data: both fits vs truth",
     ylab = "g(r)", xlab = "r", ylim = range(c(1, test_obs$iso, test_fitted_line), na.rm = TRUE))

lines(r, test_fitted_line, col = "red", lwd = 2)
lines(r, test_fitted_line2, col = "blue", lwd = 2, lty = "dashed")












