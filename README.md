# Autotune
Repository for R package of Autotune LASSO

## Installation

You can install the development version from GitHub:

```r
# Install the devtools or remotes package if you don't have it
install.packages("devtools")

devtools::install_github("Tathagata-S/Autotune")
```

👉 **Usage**:
```md
## Usage

Here’s a quick example:
```

```r
library(Autotune)
?autotune_lasso
?autotune_lasso_active
?autotune_lasso_l2

set.seed(10)
n <- 80
p <- 500
s <- 5
type <- 1
snr <- 3
betatrue <- c(rep(1,s), rep(0, p - s))
x <- scale(matrix(rnorm(n * p), ncol = p))
error.sd <- sqrt((betatrue %*% betatrue)/snr)

err <- rnorm(n, sd = error.sd)
y <- x %*% betatrue + err
y <- y - mean(y)

#Try out different functions of the package in different setups
ans <- autotune_lasso(x, y, verbose = T)
# ans <- autotune_lasso_active(x, y, verbose = T)
# ans <- autotune_lasso_l2(x, y, verbose = T)

b <- betatrue
cat("The Predictors which are actually significant:", which(b != 0), "\n")
cat("Top 10 predictors X_i's in the ranking of X_i's given by autotune:", (ans$sorted_predictors[1:10] + 1))
cat("No of significant predictors in each CD iteration when sigma_hat is allowed to vary:", ans$count_sig_beta)
cat("Sigma estimates in each CD iteration:", ans$sigma2_seq)
cat("Empirical noise variance:", noise.sd^2)
which(ans$beta != 0)
which(b != 0)

```
