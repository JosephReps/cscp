# CSCP
Fitting and simulation of LGCP's &amp; CSCP's 

## Notes

See my thesis notes here: https://josephreps.github.io/cscp/

## Visualization

See my live Cox process simulation visualizer here: https://josephreps.github.io/spatLab/

## Installation

Use

```r
devtools::install_github("https://github.com/JosephReps/cscp")
```

## Example

Simulate a CSCP:

```r
library(spatstat)

W <- owin(c(0, 10), c(0, 10))

set.seed(2)
sim <- sim_cscp(
    phi = 1, lambda = 10, s = 1,
    W = W,
    df = 1
    )
```

Plot the simulated data:

```r
par(mfrow = c(1, 2), mar = c(1.5,1.5,1.5,1.5))
plot(sim$grf[[1]], main = "Generated\n GRF")
plot(sim$intensity, main = "Corresponding\n CSCP Intensity Surface")
```

Fit a MC model:

```r
cscp_pp <- sim$pp
cscp_fit <- fit_cscp_mc(cscp_pp, bw = 0.02)
```

Plot the fitted PCF and compare to the empirical / truth:

```r
r_vals <- seq(0, 2.5, by = 0.001)

fitted_pcf_vals <- pcf_cscp_theoretical(r_vals, lambda = cscp_fit$par[["lambda"]], 
                                        phi = cscp_fit$par[["phi"]], 
                                        s = cscp_fit$par[["s"]])
theoretical_pcf_vals <- pcf_cscp_theoretical(r_vals, lambda = 100, phi = 1.5, s = 0.05)

plot(pcf(cscp_pp, correction = "isotropic", divisor = "a",
                        zerocor = "JonesFoster", bw = 0.02), 
     ylim = c(0.8, 3), legend = F, main = "Example of MC CSCP fit (lambda = 100, phi = 1.5, s = 0.05)")
lines(r_vals, theoretical_pcf_vals, col = "forestgreen", lwd = 2)
lines(r_vals, fitted_pcf_vals, col = "red", lwd = 2)
legend("topright", legend = c("Empirical", "True", "Fitted"), 
            col = c("black", "forestgreen", "red"), lwd = c(1, 2, 2))
```
