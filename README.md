# Autotune
Repository for R package of Autotune LASSO

## Installation
You can install the development version from GitHub:

```r
# Install the devtools or remotes package if you don't have it
install.packages("devtools")

# Ensure that you have the Rcpp package installed with version >=1.0.13
devtools::install_github("Tathagata-S/Autotune")
```
When installing from GitHub, in order to build the package from source, you need to have the appropriate R development tools installed ([Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows, or [these tools](https://mac.r-project.org/tools/) on Mac).


👉 **Usage**:
```md
Here’s a quick example:
```

```r
library(Autotune)
?autotune_lasso

set.seed(10)
n <- 80
p <- 400
s <- 5
snr <- 4
betatrue <- c(rep(1,s), rep(0, p - s))
x <- matrix(rnorm(n * p), ncol = p)
error.sd <- sqrt((betatrue %*% betatrue)/snr)

err <- rnorm(n, sd = error.sd)
y <- x %*% betatrue + err
y <- y - mean(y)

ans <- autotune_lasso(x, y, trace_it = T)

b <- betatrue
# The Predictors which are actually significant:
which(b != 0)
# The Predictors which had nonzero estimated coefficients:
which(ans$beta != 0)
# Top 10 predictors X_i's in the ranking of X_i's given by autotune:
ans$CD.path.details$sorted_predictors[1:10]
# No of significant predictors in each CD iteration when sigma_hat is allowed to vary:
ans$CD.path.details$count_sig_beta
# Sigma estimates in each CD iteration:
ans$CD.path.details$sigma_sq_seq
# Empirical noise variance:
var(err)
```

