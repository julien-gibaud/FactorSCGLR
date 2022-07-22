#' @title Information criterion for scglr with factors
#' @param x an object from FactorSCGLR
#' @return \item{aic}{akaike information criterion}
#' @return \item{bic}{bayesian information criterion}
#' @export
#' @examples
#' # load sample data
#' data <- genus

#'# get variable names from dataset
#'n <- names(data)
#'ny <- n[grep("^gen",n)]
#'nx1 <- n[grep("^evi",n)]
#'nx2 <- n[grep("^pluvio",n)]
#'na <- c("geology", "altitude", "forest", "lon", "lat")

#'# build multivariate formula
#'form <- multivariateFormula(Y = ny, X = list(nx1, nx2), A = na)

#'# define family
#'fam <- rep("poisson", length(ny))

#'# run function
#'H <- c(2,2)
#'J <- 2
#'met <- methodSR(l=4, s=0.5)
#'res <- FactorSCGLR(formula=form, data=data, H=H, J=J,
#'                    family=fam, method=met, offset = data$surface)
#'
#'# compute the information criterion
#'IC <- InformationCriterion(res)
#'IC$aic
#'IC$bic

InformationCriterion <- function(x){
  if (class(x) != "FactorSCGLR")
    stop("This plot function need an FactorSCGLR result")

  res <- x
  # linear predictor
  prediction <- cbind(1, res$A, res$comp) %*% res$coef

  likelihood <- 0
  # variance covariance matrix
  if(!is.null(res$B)){
    BB <- crossprod(res$B)
  } else {
    BB <- matrix(data=0, ncol=ncol(res$Y), nrow=ncol(res$Y))
  }

  # compute the log-likelihood
  for(n in 1:nrow(res$Y)){
    Sig <- BB + diag(1/res$W[n,])
    likelihood <- likelihood + Rfast::dmvnorm(x=res$Z[n,], mu=prediction[n,], sigma=Sig, logged=TRUE)
  }

  # compute the number of parameters
  npar <- 0
  npar <- npar + nrow(res$coef)*ncol(res$coef) + sum(res$family=="gaussian")
  if(!is.null(res$B)){
    npar <- npar + nrow(res$B)*ncol(res$B) - nrow(res$B)*(nrow(res$B)-1)/2
  }

  return(list(aic=-2*likelihood+2*npar,
              bic=-2*likelihood+npar*log(nrow(res$Y))))
}
