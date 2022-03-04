# weighted scale functions
# @param x vector to scale
# @param w weight
# @return scaled vector
wtScale <-function(x,w) {
  xc <- x-sum(w*x)
  v <- sum(xc^2*w)
  xcr <- xc/sqrt(v)
  return(xcr)
}

# @param K number of responses
# @return maximum number of factors
choix.num.param <- function(K){#K reponses
  return(floor((2*K+1-sqrt(8*K+1))/2))
}

# LogSumExp function
# @param m vector
# @return LogSumExp function used to compute log-likelihood
sumlogs <- function(m){
  M <- max(m)
  return(M+log(sum(exp(m-M))))
}
