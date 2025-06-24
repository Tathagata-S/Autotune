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

```r
library(Autotune)

# Run the automated parameter tuning
results <- autotune(data = my_data, ...)
head(results)


