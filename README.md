# FactorsSCGLR
Supervised Component-based Generalized Linear Regression for factor models

## Installation

``` r
# Install development version from GitHub
remotes::install_github("julien-gibaud/FactorsSCGLR")
```

## Example
``` r
# load sample data
data <- genus
data <- as.data.frame(apply(data, 2, as.numeric ))

# get variable names from dataset
n <- names(data)
ny <- n[grep("^gen",n)]
nx1 <- n[grep("^evi",n)]
nx2 <- n[grep("^pluvio",n)]
na <- c("geology", "altitude", "forest", "lon", "lat")

# build multivariate formula
form <- multivariateFormula(Y = ny, X = list(nx1, nx2), A = na)

# define family
fam <- rep("poisson", length(ny))

# run function
H <- c(2,2)
J <- 2
met <- methodSR(l=4, s=0.5)
res <- FactorsSCGLR(formula=form, data=data, H=H, J=J,
                    family=fam, method=met, offset = data$surface)
                  
# print results
res$U
res$comp
res$B

# plot the results
plot_comp(x=res, thresold=0.5, theme=1, plan=c(1,2))
plot_comp(x=res, thresold=0.5, theme=2, plan=c(1,2))
```
