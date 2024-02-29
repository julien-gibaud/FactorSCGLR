# FactorSCGLR
Supervised Component-based Generalized Linear Regression for factor models

## Installation

``` r
# Install development version from GitHub
remotes::install_github("julien-gibaud/FactorSCGLR")
```

## Example
``` r
library(FactorSCGLR)

# load sample data
data <- genus

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
res <- FactorSCGLR(formula=form, data=data, H=H, J=J,
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

# Detect the clusters 
CD <- ClusterDetection(mat=res$B)
#the identified clusters
CD$cluster
#the output of the multidimensional scaling
CD$mds
```

## Simulations from the paper
### Simulation with a mix of distribution
``` r
library(FactorSCGLR)
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
  Y <- cbind(Y, rnorm(n=N, mean=mu, sd=1))
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
res <- FactorSCGLR( formula=form,
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
CD <- ClusterDetection(mat=res$B)
#the identified clusters
CD$cluster
#the output of the multidimensional scaling
CD$mds

```

### Packages comparison
#### Gaussian distribution
```r
library(FactorSCGLR)
library(gllvm)
library(mvtnorm)
library(fossil)
library(vegan)


#**********#
# settings #
#**********#
N <- 100     # observations
K <- 10      # responses
J <- 2       # factors
variance_B  <- 0.1 # variance within the clusters

#*****************************#
# create the latent variables #
#*****************************#
psi1.sim <- rnorm(n = N)
psi2.sim <- rnorm(n = N)
G.sim <- rmvnorm(n=N, mean = rep(0, J), sigma = diag(J))

#*********************************#
# create the additional covariate #
#*********************************#
A <- as.factor(sample(x = 1:8, size = N, replace = TRUE))
while (length(levels(A)) < 8) as.factor(sample(x = 1:8, size = N, replace = TRUE))
A_ind <- model.matrix(~-1+A)

#**********************************#
# create the regression parameters #
#**********************************#
gamma1.sim <- runif(n = K, min = -4, max = 4)
gamma2.sim <- runif(n = K, min = -2, max = 2)
delta.sim <- matrix(data = runif(n = 8*K, min = -0.5, max = 0.5), nrow = 8, ncol = K)

#*****************************#
# create the factors loadings #
#*****************************#
mean1 <- c(0,2)
mean2 <- c(1.5,0)
rot <- rep(-1, K)
rot[(1:K)%%2==0] <- 1
B.sim <- t(diag(rot)%*%rbind(rmvnorm(n=0.4*K, mean = mean1, sigma = variance_B*diag(J)),
                             rmvnorm(n=0.6*K, mean = mean2, sigma = variance_B*diag(J))))

#*********************************#
# simulate the response variables #
#*********************************#
Y <- cbind()
for(i in 1:K){ 
  eta <- psi1.sim*gamma1.sim[i]+psi2.sim*gamma2.sim[i]+A_ind%*%delta.sim[,i]+G.sim%*%B.sim[,i]
  mu <- eta
  Y <- cbind(Y, rnorm(n=N, mean=mu, sd=1))
}

#************************************#
# simulate the explanatory variables #
#************************************#
X <- cbind()
for(i in 1:5)  X <- cbind(X, psi1.sim + rnorm(n = N, sd = sqrt(0.1)))
for(i in 1:5)  X <- cbind(X, psi2.sim + rnorm(n = N, sd = sqrt(0.1)))
for(i in 1:10) X <- cbind(X, rnorm(n = N))

#**************#
# run function #
#**************#
# build data
data <- data.frame(cbind(Y,X))
data$A <- A
# build multivariate formula
ny <- paste("Y", 1:K, sep = "")
nx <- paste("X", 1:ncol(X), sep = "")
na <- "A"
colnames(data) <- c(ny,nx,na)
form <- multivariateFormula(ny,nx,na, additional = TRUE)
# define family
fam <- rep('gaussian', K)
# define method
met <- methodSR(l=4, s=0.5)
# run FactorsSCGLR
H <- 2
start.time.scglr <- Sys.time()
res.scglr <- FactorSCGLR( formula=form,
                          data=data,
                          J=J,
                          H=H,
                          method=met,
                          family = fam)
end.time.scglr <- Sys.time()
time.scglr <- as.numeric(difftime(end.time.scglr,
                                  start.time.scglr, 
                                  units = "secs"))

# Detect the clusters 
cluster.scglr <- ClusterDetection(mat=res.scglr$B)
# Compute RI
RI.scglr <- rand.index(cluster.scglr$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Compute ARI
ARI.scglr <- adj.rand.index(cluster.scglr$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Procrustes errors
proc.B.scglr <- procrustes(t(B.sim), t(res.scglr$B))$ss
proc.G.scglr <- procrustes(G.sim, res.scglr$G)$ss
# Latent dimension recovery
cor1.scglr <- max(cor(psi1.sim, res.scglr$comp)^2)
cor2.scglr <- max(cor(psi2.sim, res.scglr$comp)^2)

#***********#
# run gllvm #
#***********#
# build data
explanatory <- cbind(X,A_ind[,-1])
colnames(explanatory) <- c(nx, paste("A", 2:8, sep = ""))
# run gllvm-VA
start.time.gllvm.va <- Sys.time()
res.gllvm.va <- gllvm(Y,
                      X = explanatory, 
                      formula = ~ A2+A3+A4+A5+A6+A7+A8,
                      family = gaussian(),
                      num.lv = 2,
                      num.RR = 2,
                      lv.formula = ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20,
                      method = "VA",
                      sd.errors = FALSE)
end.time.gllvm.va <- Sys.time()
time.gllvm.va <- as.numeric(difftime(end.time.gllvm.va,
                                     start.time.gllvm.va,
                                     units = "secs"))

# Detect the clusters 
B.gllvm.va <- t(res.gllvm.va$params$theta[,3:4]%*%diag(res.gllvm.va$params$sigma.lv,2))
cluster.gllvm.va <- ClusterDetection(B.gllvm.va)
# Compute RI
RI.gllvm.va <- rand.index(cluster.gllvm.va$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Compute ARI
ARI.gllvm.va <- adj.rand.index(cluster.gllvm.va$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Procrustes errors
proc.B.gllvm.va <- procrustes(t(B.sim), t(B.gllvm.va))$ss
factors.gllvm.va <- getLV(res.gllvm.va, type = 'residual')
proc.G.gllvm.va <- procrustes(G.sim, factors.gllvm.va)$ss
# Latent dimension recovery
comps.gllvm.va <- getLV(res.gllvm.va, type = 'marginal')
cor1.gllvm.va <- max(cor(psi1.sim, comps.gllvm.va)^2)
cor2.gllvm.va <- max(cor(psi2.sim, comps.gllvm.va)^2)


# run gllvm-LA
start.time.gllvm.la <- Sys.time()
res.gllvm.la <- gllvm(Y,
                      X = explanatory,
                      formula = ~ A2+A3+A4+A5+A6+A7+A8,
                      family = gaussian(),
                      num.lv = 2,
                      num.RR = 2,
                      lv.formula = ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20,
                      method = "LA",
                      sd.errors = FALSE)
end.time.gllvm.la <- Sys.time()
time.gllvm.la <- as.numeric(difftime(end.time.gllvm.la,
                                     start.time.gllvm.la,
                                     units = "secs"))

# Detect the clusters 
B.gllvm.la <- t(res.gllvm.la$params$theta[,3:4]%*%diag(res.gllvm.la$params$sigma.lv,2))
cluster.gllvm.la <- ClusterDetection(B.gllvm.la)
# Compute RI
RI.gllvm.la <- rand.index(cluster.gllvm.la$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Compute ARI
ARI.gllvm.la <- adj.rand.index(cluster.gllvm.la$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Procrustes errors
proc.B.gllvm.la <- procrustes(t(B.sim), t(B.gllvm.la))$ss
factors.gllvm.la <- getLV(res.gllvm.la, type = 'residual')
proc.G.gllvm.la <- procrustes(G.sim, factors.gllvm.la)$ss
# Latent dimension recovery
comps.gllvm.la <- getLV(res.gllvm.la, type = 'marginal')
cor1.gllvm.la <- max(cor(psi1.sim, comps.gllvm.la)^2)
cor2.gllvm.la <- max(cor(psi2.sim, comps.gllvm.la)^2)
```
#### Poisson distribution
```r
library(FactorSCGLR)
library(gllvm)
library(mvtnorm)
library(fossil)
library(vegan)


#**********#
# settings #
#**********#
N <- 100     # observations
K <- 10      # responses
J <- 2       # factors
variance_B  <- 0.1 # variance within the clusters

#*****************************#
# create the latent variables #
#*****************************#
psi1.sim <- rnorm(n = N)
psi2.sim <- rnorm(n = N)
G.sim <- rmvnorm(n=N, mean = rep(0, J), sigma = diag(J))

#*********************************#
# create the additional covariate #
#*********************************#
A <- as.factor(sample(x = 1:8, size = N, replace = TRUE))
while (length(levels(A)) < 8) as.factor(sample(x = 1:8, size = N, replace = TRUE))
A_ind <- model.matrix(~-1+A)

#**********************************#
# create the regression parameters #
#**********************************#
gamma1.sim <- runif(n = K, min = -4, max = 4)
gamma2.sim <- runif(n = K, min = -2, max = 2)
delta.sim <- matrix(data = runif(n = 8*K, min = -0.5, max = 0.5), nrow = 8, ncol = K)

#*****************************#
# create the factors loadings #
#*****************************#
mean1 <- c(0,2)
mean2 <- c(1.5,0)
rot <- rep(-1, K)
rot[(1:K)%%2==0] <- 1
B.sim <- t(diag(rot)%*%rbind(rmvnorm(n=0.4*K, mean = mean1, sigma = variance_B*diag(J)),
                             rmvnorm(n=0.6*K, mean = mean2, sigma = variance_B*diag(J))))

#*********************************#
# simulate the response variables #
#*********************************#
Y <- cbind()
for(i in 1:K){
  eta <- 0.75*(psi1.sim*gamma1.sim[i]+psi2.sim*gamma2.sim[i]+A_ind%*%delta.sim[,i]+G.sim%*%B.sim[,i])
  mu <- exp(eta)
  Y <- cbind(Y, rpois(n = N, lambda = mu))
}

#************************************#
# simulate the explanatory variables #
#************************************#
X <- cbind()
for(i in 1:5)  X <- cbind(X, psi1.sim + rnorm(n = N, sd = sqrt(0.1)))
for(i in 1:5)  X <- cbind(X, psi2.sim + rnorm(n = N, sd = sqrt(0.1)))
for(i in 1:10) X <- cbind(X, rnorm(n = N))

#**************#
# run function #
#**************#
# build data
data <- data.frame(cbind(Y,X))
data$A <- A
# build multivariate formula
ny <- paste("Y", 1:K, sep = "")
nx <- paste("X", 1:ncol(X), sep = "")
na <- "A"
colnames(data) <- c(ny,nx,na)
form <- multivariateFormula(ny,nx,na, additional = TRUE)
# define family
fam <- rep('poisson', K)
# define method
met <- methodSR(l=4, s=0.5)
# run FactorsSCGLR
H <- 2
start.time.scglr <- Sys.time()
res.scglr <- FactorSCGLR( formula=form,
                          data=data,
                          J=J,
                          H=H,
                          method=met,
                          family = fam)
end.time.scglr <- Sys.time()
time.scglr <- as.numeric(difftime(end.time.scglr,
                                  start.time.scglr, 
                                  units = "secs"))

# Detect the clusters 
cluster.scglr <- ClusterDetection(mat=res.scglr$B)
# Compute RI
RI.scglr <- rand.index(cluster.scglr$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Compute ARI
ARI.scglr <- adj.rand.index(cluster.scglr$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Procrustes errors
proc.B.scglr <- procrustes(t(B.sim), t(res.scglr$B))$ss
proc.G.scglr <- procrustes(G.sim, res.scglr$G)$ss
# Latent dimension recovery
cor1.scglr <- max(cor(psi1.sim, res.scglr$comp)^2)
cor2.scglr <- max(cor(psi2.sim, res.scglr$comp)^2)

#***********#
# run gllvm #
#***********#
# build data
explanatory <- cbind(X,A_ind[,-1])
colnames(explanatory) <- c(nx, paste("A", 2:8, sep = ""))
# run gllvm-VA
start.time.gllvm.va <- Sys.time()
res.gllvm.va <- gllvm(Y,
                      X = explanatory, 
                      formula = ~ A2+A3+A4+A5+A6+A7+A8,
                      family = poisson(),
                      num.lv = 2,
                      num.RR = 2,
                      lv.formula = ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20,
                      method = "VA",
                      sd.errors = FALSE)
end.time.gllvm.va <- Sys.time()
time.gllvm.va <- as.numeric(difftime(end.time.gllvm.va,
                                     start.time.gllvm.va,
                                     units = "secs"))

# Detect the clusters 
B.gllvm.va <- t(res.gllvm.va$params$theta[,3:4]%*%diag(res.gllvm.va$params$sigma.lv,2))
cluster.gllvm.va <- ClusterDetection(B.gllvm.va)
# Compute RI
RI.gllvm.va <- rand.index(cluster.gllvm.va$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Compute ARI
ARI.gllvm.va <- adj.rand.index(cluster.gllvm.va$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Procrustes errors
proc.B.gllvm.va <- procrustes(t(B.sim), t(B.gllvm.va))$ss
factors.gllvm.va <- getLV(res.gllvm.va, type = 'residual')
proc.G.gllvm.va <- procrustes(G.sim, factors.gllvm.va)$ss
# Latent dimension recovery
comps.gllvm.va <- getLV(res.gllvm.va, type = 'marginal')
cor1.gllvm.va <- max(cor(psi1.sim, comps.gllvm.va)^2)
cor2.gllvm.va <- max(cor(psi2.sim, comps.gllvm.va)^2)


# run gllvm-LA
start.time.gllvm.la <- Sys.time()
res.gllvm.la <- gllvm(Y,
                      X = explanatory,
                      formula = ~ A2+A3+A4+A5+A6+A7+A8,
                      family = poisson(),
                      num.lv = 2,
                      num.RR = 2,
                      lv.formula = ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20,
                      method = "LA",
                      sd.errors = FALSE)
end.time.gllvm.la <- Sys.time()
time.gllvm.la <- as.numeric(difftime(end.time.gllvm.la,
                                     start.time.gllvm.la,
                                     units = "secs"))

# Detect the clusters 
B.gllvm.la <- t(res.gllvm.la$params$theta[,3:4]%*%diag(res.gllvm.la$params$sigma.lv,2))
cluster.gllvm.la <- ClusterDetection(B.gllvm.la)
# Compute RI
RI.gllvm.la <- rand.index(cluster.gllvm.la$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Compute ARI
ARI.gllvm.la <- adj.rand.index(cluster.gllvm.la$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Procrustes errors
proc.B.gllvm.la <- procrustes(t(B.sim), t(B.gllvm.la))$ss
factors.gllvm.la <- getLV(res.gllvm.la, type = 'residual')
proc.G.gllvm.la <- procrustes(G.sim, factors.gllvm.la)$ss
# Latent dimension recovery
comps.gllvm.la <- getLV(res.gllvm.la, type = 'marginal')
cor1.gllvm.la <- max(cor(psi1.sim, comps.gllvm.la)^2)
cor2.gllvm.la <- max(cor(psi2.sim, comps.gllvm.la)^2)
```
#### Bernoulli distribution
```r
library(FactorSCGLR)
library(gllvm)
library(mvtnorm)
library(fossil)
library(vegan)


#**********#
# settings #
#**********#
N <- 100     # observations
K <- 10      # responses
J <- 2       # factors
variance_B  <- 0.1 # variance within the clusters
set.seed(42)

#*****************************#
# create the latent variables #
#*****************************#
psi1.sim <- rnorm(n = N)
psi2.sim <- rnorm(n = N)
G.sim <- rmvnorm(n=N, mean = rep(0, J), sigma = diag(J))

#*********************************#
# create the additional covariate #
#*********************************#
A <- as.factor(sample(x = 1:8, size = N, replace = TRUE))
while (length(levels(A)) < 8) as.factor(sample(x = 1:8, size = N, replace = TRUE))
A_ind <- model.matrix(~-1+A)

#**********************************#
# create the regression parameters #
#**********************************#
gamma1.sim <- runif(n = K, min = -4, max = 4)
gamma2.sim <- runif(n = K, min = -2, max = 2)
delta.sim <- matrix(data = runif(n = 8*K, min = -0.5, max = 0.5), nrow = 8, ncol = K)

#*****************************#
# create the factors loadings #
#*****************************#
mean1 <- c(0,2)
mean2 <- c(1.5,0)
rot <- rep(-1, K)
rot[(1:K)%%2==0] <- 1
B.sim <- t(diag(rot)%*%rbind(rmvnorm(n=0.4*K, mean = mean1, sigma = variance_B*diag(J)),
                             rmvnorm(n=0.6*K, mean = mean2, sigma = variance_B*diag(J))))

#*********************************#
# simulate the response variables #
#*********************************#
Y <- cbind()
for(i in 1:K){
  eta <- psi1.sim*gamma1.sim[i]+psi2.sim*gamma2.sim[i]+A_ind%*%delta.sim[,i]+G.sim%*%B.sim[,i]
  mu <- exp(eta)/(1+exp(eta))
  Y <- cbind(Y, rbinom(n=N, size=1, prob=mu))
}

#************************************#
# simulate the explanatory variables #
#************************************#
X <- cbind()
for(i in 1:5)  X <- cbind(X, psi1.sim + rnorm(n = N, sd = sqrt(0.1)))
for(i in 1:5)  X <- cbind(X, psi2.sim + rnorm(n = N, sd = sqrt(0.1)))
for(i in 1:10) X <- cbind(X, rnorm(n = N))

#**************#
# run function #
#**************#
# build data
data <- data.frame(cbind(Y,X))
data$A <- A
# build multivariate formula
ny <- paste("Y", 1:K, sep = "")
nx <- paste("X", 1:ncol(X), sep = "")
na <- "A"
colnames(data) <- c(ny,nx,na)
form <- multivariateFormula(ny,nx,na, additional = TRUE)
# define family
fam <- rep('bernoulli', K)
# define method
met <- methodSR(l=4, s=0.5)
# run FactorsSCGLR
H <- 2
start.time.scglr <- Sys.time()
res.scglr <- FactorSCGLR( formula=form,
                          data=data,
                          J=J,
                          H=H,
                          method=met,
                          family = fam)
end.time.scglr <- Sys.time()
time.scglr <- as.numeric(difftime(end.time.scglr,
                                  start.time.scglr, 
                                  units = "secs"))

# Detect the clusters 
cluster.scglr <- ClusterDetection(mat=res.scglr$B)
# Compute RI
RI.scglr <- rand.index(cluster.scglr$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Compute ARI
ARI.scglr <- adj.rand.index(cluster.scglr$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Procrustes errors
proc.B.scglr <- procrustes(t(B.sim), t(res.scglr$B))$ss
proc.G.scglr <- procrustes(G.sim, res.scglr$G)$ss
# Latent dimension recovery
cor1.scglr <- max(cor(psi1.sim, res.scglr$comp)^2)
cor2.scglr <- max(cor(psi2.sim, res.scglr$comp)^2)

#***********#
# run gllvm #
#***********#
# build data
explanatory <- cbind(scale(X),A_ind[,-1])
colnames(explanatory) <- c(nx, paste("A", 2:8, sep = ""))
# run gllvm-EVA
start.time.gllvm.eva <- Sys.time()
res.gllvm.eva <- gllvm(Y,
                       X = explanatory,
                       formula = ~ A2+A3+A4+A5+A6+A7+A8,
                       family = binomial(link = "logit"),
                       num.lv = 2,
                       num.RR = 2,
                       lv.formula = ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20,
                       method = "EVA",
                       sd.errors = FALSE)
end.time.gllvm.eva <- Sys.time()
time.gllvm.eva <- as.numeric(difftime(end.time.gllvm.eva,
                                      start.time.gllvm.eva,
                                      units = "secs"))

# Detect the clusters 
B.gllvm.eva <- t(res.gllvm.eva$params$theta[,3:4]%*%diag(res.gllvm.eva$params$sigma.lv,2))
cluster.gllvm.eva <- ClusterDetection(B.gllvm.eva)
# Compute RI
RI.gllvm.eva <- rand.index(cluster.gllvm.eva$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Compute ARI
ARI.gllvm.eva <- adj.rand.index(cluster.gllvm.eva$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Procrustes errors
proc.B.gllvm.eva <- procrustes(t(B.sim), t(B.gllvm.eva))$ss
factors.gllvm.eva <- getLV(res.gllvm.eva, type = 'residual')
proc.G.gllvm.eva <- procrustes(G.sim, factors.gllvm.eva)$ss
# Latent dimension recovery
comps.gllvm.eva <- getLV(res.gllvm.eva, type = 'marginal')
cor1.gllvm.eva <- max(cor(psi1.sim, comps.gllvm.eva)^2)
cor2.gllvm.eva <- max(cor(psi2.sim, comps.gllvm.eva)^2)


# run gllvm-LA
start.time.gllvm.la <- Sys.time()
res.gllvm.la <- gllvm(Y,
                      X = explanatory,
                      formula = ~ A2+A3+A4+A5+A6+A7+A8,
                      family = binomial(link = "logit"),
                      num.lv = 2,
                      num.RR = 2,
                      lv.formula = ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20,
                      method = "LA",
                      sd.errors = FALSE)
end.time.gllvm.la <- Sys.time()
time.gllvm.la <- as.numeric(difftime(end.time.gllvm.la,
                                     start.time.gllvm.la,
                                     units = "secs"))

# Detect the clusters 
B.gllvm.la <- t(res.gllvm.la$params$theta[,3:4]%*%diag(res.gllvm.la$params$sigma.lv,2))
cluster.gllvm.la <- ClusterDetection(B.gllvm.la)
# Compute RI
RI.gllvm.la <- rand.index(cluster.gllvm.la$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Compute ARI
ARI.gllvm.la <- adj.rand.index(cluster.gllvm.la$cluster, c(rep(1,0.4*K),rep(2,0.6*K)))
# Procrustes errors
proc.B.gllvm.la <- procrustes(t(B.sim), t(B.gllvm.la))$ss
factors.gllvm.la <- getLV(res.gllvm.la, type = 'residual')
proc.G.gllvm.la <- procrustes(G.sim, factors.gllvm.la)$ss
# Latent dimension recovery
comps.gllvm.la <- getLV(res.gllvm.la, type = 'marginal')
cor1.gllvm.la <- max(cor(psi1.sim, comps.gllvm.la)^2)
cor2.gllvm.la <- max(cor(psi2.sim, comps.gllvm.la)^2)

```

### Application to a real dataset
```r
library(FactorSCGLR)
library(SCGLR)

# load data available at https://doi.org/10.15454/AJZUQN
load('Copie de data_Duflot et al_AGEE.RData')
data <- x
data$year <- as.factor(data$year)

# get variable names from dataset
n <- names(data)

# agrobiodiversity (responses)
ny <- c("carab.richn.tot", "carab.abund.tot", "carab.shannon.tot", "c.ax1.coa",
        "c.ax2.coa", "c.ax3.coa", "plants.richn", "plants.abund",
        "plants.shannon", "p.ax1.coa", "p.ax2.coa", "p.ax3.coa")

# pest control (first theme)
nx1 <- c("aphid.low.tot", "aphid.high.tot", "seeds.tot", "eggs.tot")
# farming intensity (second theme)
nx2 <- c("TFI.f", "TFI.h", "TFI.total", "qtyN.kg", "cum.till.depth", "nb.op")
# land. hetero. SNC (third theme)
nx3 <- c("X.SNC", "X.Pgrass", "X.Wooded", "MPSSNH", "edgesSNHcrop", "edgesSNH")
# land. hetero. crop mosaic (fourth theme)
nx4 <- c("X.W.cereral", "X.otherWcrop", "X.S.crop", "SHDIcrop", "edgedensitycrop")
# all explanatory variables
nx <- c(nx1, nx2, nx3, nx4)

# additional covariate
na <- c("year")

# define family
fam <- c('poisson', 'poisson', 'gaussian', 'gaussian', 'gaussian', 'gaussian',
         'poisson', 'gaussian', 'gaussian', 'gaussian', 'gaussian', 'gaussian')

#*********************************************#
# calibration of the hyper-parameters s and l # 
#*********************************************#

# build multivariate formula for SCGLR
form <- multivariateFormula(Y = ny, X = nx, A = na)

# we print the combination minimizing the cross-validation
S <- c(0.1, 0.3, 0.5)
L <- c(1, 2, 3, 4, 7, 10)
folds <- c(rep(1, 10), rep(2, 10), rep(3, 10), rep(4, 10), rep(5, 14))
min_cv <- Inf
for(s in S){
  for(l in L){
    genus.cv <- scglrCrossVal(formula=form, data=data, family=fam, K=3,
                              method=methodSR(l=l,s=s),
                              folds = folds)
    mean.crit <- colMeans(log(genus.cv))
    cv <- mean(mean.crit)
    if(cv < min_cv){
      min_cv <- cv
      print(paste("s=", s, "l=", l, "CV=", cv))
    }
  }
}
# we keep s = 0.5 and l = 1
 
#********************************************************#
# calibration of the hyper-parameters H1, H2, H3, and H4 # 
#********************************************************#

# build multivariate formula for F-SCGLR
form <- multivariateFormula(Y = ny, X = list(nx1, nx2, nx3, nx4), A = na)

# we print the combination minimizing the BIC
min_bic <- Inf
H <- 3
for(h1 in 0:H){
  for(h2 in 0:H){
    for(h3 in 0:H){
      for(h4 in 0:H){
        res <-  FactorSCGLR(formula=form, data=data, J=0,
                            H=c(h1,h2,h3,h4),
                            method=methodSR(l=1,s=0.5),
                            family = fam)
        crit <- InformationCriterion(x=res)
        if(crit$bic < min_bic){
          min_bic <- crit$bic
          print(paste("H1=", h1, "H2=", h2, "H3=", h3, "H4=", h4,
                      "bic=", round(crit$bic, digits = 0)))
        }
      }
    }
  }
}
# we keep (H1, H2, H3, H4) = (0, 3, 0, 0)

#**************************************#
# calibration of the hyper-parameter J # 
#**************************************#

# we print the value of J minimizing the BIC
min_bic <- Inf
for(j in 0:5){
  res <-  FactorSCGLR(formula=form, data=data, J=j,
                      H=c(0,3,0,0),
                      method=methodSR(l=1,s=0.5),
                      family = fam)
  crit <- InformationCriterion(x=res)
  if(crit$bic < min_bic){
    min_bic <- crit$bic
    print(paste("J=", j, "bic=", round(crit$bic, digits = 0)))
  }
}
# we keep J = 3

#*********************#
# run the final model #
#*********************#

# run
res <-  FactorSCGLR(formula=form, data=data,
                    J=3, H=c(0,3,0,0),
                    method=methodSR(l=1,s=0.5),
                    family = fam)

# the correlation plots
plot1 <- plot_comp(x=res, thresold = 0.75, theme=2, plan = c(1,2))
plot2 <- plot_comp(x=res, thresold = 0.75, theme=2, plan = c(1,3))
plot3 <- plot_comp(x=res, thresold = 0.75, theme=2, plan = c(2,3))
# the supervised components
res$comp
# the factor loadings
res$B
# the factors
res$G
# the residual variances of the Gaussian responses
res$sigma2

# Detect the clusters 
CD <- ClusterDetection(mat=res$B)
#the identified clusters
CD$cluster
#the output of the multidimensional scaling
CD$mds
```
