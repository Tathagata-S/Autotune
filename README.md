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

👉 **Usage**:
```md
Here’s a quick example:
```

```r
library(Autotune)
?autotune_lasso

set.seed(1234)
n <- 80
p <- 500
s <- 5
snr <- 3
betatrue <- c(rep(1,s), rep(0, p - s))
x <- matrix(rnorm(n * p), ncol = p)
error.sd <- sqrt((betatrue %*% betatrue)/snr)

err <- rnorm(n, sd = error.sd)
y <- x %*% betatrue + err
y <- y - mean(y)

ans <- autotune_lasso(x, y, verbose = T)

b <- betatrue
# The Predictors which are actually significant:
which(b != 0)
# The Predictors which had nonzero estmated coefficients:
which(ans$beta != 0)
# Top 10 predictors X_i's in the ranking of X_i's given by autotune:
ans$sorted_predictors[1:10] + 1
# No of significant predictors in each CD iteration when sigma_hat is allowed to vary:
ans$count_sig_beta
# Sigma estimates in each CD iteration:
ans$sigma2_seq
# Empirical noise variance:
var(err)
```
