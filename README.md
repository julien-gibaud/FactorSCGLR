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
na <- c("geology", "forest")

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
res$coef

# plot the results
plot_comp(x=res, thresold=0.5, theme=1, plan=c(1,2))
plot_comp(x=res, thresold=0.5, theme=2, plan=c(1,2))

# compute the information criterion
IC <- InformationCriterion(res)
IC$aic
IC$bic
```

## Simulations from the paper
### Simulation with a mix of distribution
``` r
library(FactorsSCGLR)
library(mvtnorm)

#**********#
# settings #
#**********#
N <- 100     # observations
K <- 50      # responses
J <- 3       # factors
variance_B  <- 0.1 # variance within the clusters

#*****************************#
# create the latent variables #
#*****************************#
psi1.sim <- rnorm(n = N)
psi2.sim <- rnorm(n = N)
psi3.sim <- rnorm(n = N)
psi4.sim <- rnorm(n = N)
psi5.sim <- rnorm(n = N)
G.sim <- rmvnorm(n=N, mean = rep(0, J), sigma = diag(J))

#**********************************#
# create the regression parameters #
#**********************************#
gamma.sim <- runif(n = K, min = -4, max = 4)
gamma.sim1 <- runif(n = K, min = -2, max = 2)
gamma.sim2 <- runif(n = K, min = -0.5, max = 0.5)

#********************************************#
# create the variance for Gaussian responses #
#********************************************#
sigma2.sim <- runif(n = K, min = 0.1, max = 0.2)

#*****************************#
# create the factors loadings #
#*****************************#
mean1 <- c(2,0,0)
mean2 <- c(0,-1,0)
mean3 <- c(0,0,1.5)
B.sim <- t(rbind(rmvnorm(n=5, mean = mean1, sigma = variance_B*diag(J)),
                 rmvnorm(n=5, mean = -mean1, sigma = variance_B*diag(J)),
                 rmvnorm(n=10, mean = mean2, sigma = variance_B*diag(J)),
                 rmvnorm(n=15, mean = -mean2, sigma = variance_B*diag(J)),
                 rmvnorm(n=15, mean = mean3, sigma = variance_B*diag(J))))

#*********************************#
# simulate the response variables #
#*********************************#
Y <- cbind()
for(i in 1:20){ # gaussien
  eta <- psi1.sim*gamma.sim[i]+psi2.sim*gamma.sim1[i]+G.sim%*%B.sim[,i]
  mu <- eta
  Y <- cbind(Y, rnorm(n=N, mean=mu, sd=sqrt(sigma2.sim[i])))
}
for(i in 21:40){ # poisson
  eta <- 0.5*(psi4.sim*gamma.sim[i]+psi3.sim*gamma.sim1[i])+G.sim%*%B.sim[,i]
  mu <- exp(eta)
  Y <- cbind(Y, rpois(n=N, lambda=mu))
}
for(i in 41:50){ # bernoulli
  eta <- psi3.sim*gamma.sim1[i]+psi2.sim*gamma.sim2[i]+G.sim%*%B.sim[,i]
  mu <- exp(eta)/(1+exp(eta))
  Y <- cbind(Y, rbinom(n=N, size=1, prob=mu))
}

#************************************#
# simulate the explanatory variables #
#************************************#
X <- cbind()
# theme 1
for(i in 1:90)  X <- cbind(X, psi1.sim + rnorm(n = N, sd = sqrt(0.1)))
for(i in 91:150)  X <- cbind(X, psi2.sim + rnorm(n = N, sd = sqrt(0.1)))
for(i in 151:200)  X <- cbind(X, rnorm(n = N, sd = 1))
# theme 2
for(i in 1:100)  X <- cbind(X, psi3.sim + rnorm(n = N, sd = sqrt(0.1)))
for(i in 101:180)  X <- cbind(X, psi4.sim + rnorm(n = N, sd = sqrt(0.1)))
for(i in 181:240)  X <- cbind(X, psi5.sim + rnorm(n = N, sd = sqrt(0.1)))
for(i in 241:300)  X <- cbind(X, rnorm(n = N, sd = 1))

#**************#
# run function #
#**************#
# build data
data <- data.frame(cbind(Y,X))
# build multivariate formula
ny <- paste("Y", 1:50, sep = "")
nx1 <- paste("X", 1:200, sep = "") #theme 1
nx2 <- paste("X", 201:500, sep = "") #theme 2
colnames(data) <- c(ny,nx1,nx2)
form <- multivariateFormula(ny,nx1,nx2, additional = F)
# define family
fam <- c(rep('gaussian', 20), 
         rep('poisson', 20), 
         rep('bernoulli', 10))
# define method
met <- methodSR(l=4,s=0.3, maxiter = 50)
# define crit
crit <- list(tol = 1e-6, maxit = 100)
# run
H <- c(2,2)
res <- FactorsSCGLR(formula=form,
                    data=data, 
                    J=J,
                    H=H,
                    method=met, 
                    family = fam,
                    crit = crit)

#**************#
# show results #
#**************#
# the correlation plots
plot1 <- plot_comp(x=res, thresold = 0.5, theme=1, plan = c(1,2))
plot2 <- plot_comp(x=res, thresold = 0.5, them=2, plan = c(1,2))
# the supervised components
res$comp
# the factors loading
res$B
# the factors
res$G
# the residual variances of the gaussian responses
res$sigma2

#***********************************#
# compute the information criterion #
#***********************************#
IC <- InformationCriterion(x=res)
IC$aic
IC$bic

#*********************#
# Detect the clusters #
#*********************#
cluster <- ClusterDetection(res$B)
```
