#' @title Function that fits the scglr model with factors
#'
#' @param formula an object of class \code{MultivariateFormula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data a data frame to be modeled.
#' @param family a vector of character of the same length as the number of dependent variables:
#' "bernoulli", "binomial", "poisson" or "gaussian" is allowed.
#' @param H a vector of integer representing the number of components per theme.
#' @param J the number of factors
#' @param size describes the number of trials for the binomial dependent variables.
#' @param offset used for the poisson dependent variables.
#' @param crit a list of two elements : maxit and tol, describing respectively the maximum number of iterations and
#' the tolerance convergence criterion for the Expectation Maximization algorithm. Default is set to 50 and 10e-6 respectively.
#' @param method structural relevance criterion.
#' @examples
#' # load sample data
#' data <- genus
#' data <- as.data.frame(apply(data, 2, as.numeric ))

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
#'res <- FactorsSCGLR(formula=form, data=data, H=H, J=J,
#'                    family=fam, method=met, offset = data$surface)
#'
#' @return \item{U}{the set of loading vectors}
#' @return \item{comp}{the set of components}
#' @return \item{eta}{the fitted linear predictor}
#' @return \item{coef}{contains regression parameters: intercept, gamma and delta}
#' @return \item{B}{the factors loading}
#' @return \item{G}{the set of factors}
#' @return \item{sigma2}{the residual variances for gaussian responses}
#' @return \item{offset}{the offset used for the poisson dependent variables}
#' @return \item{WorkingVariables}{the set of working variables}
#' @return \item{Y}{the response matrix}
#' @return \item{Theme}{the list of standardized explanatory themes}
#' @return \item{A}{the supplementary explanatory variables}
#' @return \item{W}{the matrix of GLM weights}
#' @return \item{family}{the distributions}
#' @export
#'


FactorsSCGLR <- function(formula,
                         data,
                         family,
                         H=c(1,1),
                         J=1,
                         size=NULL,
                         offset=NULL,
                         crit=list(),
                         method=methodSR()){

  if(!inherits(formula,"MultivariateFormula"))
    formula <- multivariateFormula(formula,data=data)

  additional <- formula$additional

  # check data
  if(!inherits(data, "data.frame"))
    stop("data must be compatible with a data.frame!")
  data_size <- nrow(data)

  # check crit
  crit <- do.call("critConvergence", crit)

  # extract left-hand side (Y)
  # Y is a formula of the form ~...
  if(length(formula)[1] != 1)
    stop("Left hand side part of formula (Y) must have ONE part!")
  terms_Y <- stats::terms(formula, lhs=1, rhs=0)
  Y_vars <- all.vars(terms_Y)

  # check family
  if(!is.character(family)||any(!(family%in%c("gaussian","poisson","bernoulli","binomial"))))
    stop("Expecting character vector of gaussian, poisson, bernoulli or binomial for family")
  if(!(length(family)%in%c(1,length(Y_vars))))
    stop("Length of family must be equal to one or number of Y variables")
  if(length(family)==1)
    family <- rep(family,length(Y_vars))

  # check and preprocess size parameter
  binomial_count <- sum(family=="binomial")
  if(!is.null(size)&&(binomial_count>0)) {
    if(is.vector(size)) {
      if(sum(family=="binomial")>1)
        message("Assuming that each binomial variable has same size!")
      size <- matrix(size, data_size, binomial_count)
    } else {
      size <- as.matrix(size)
    }
  }

  # check and preprocess offset parameter
  poisson_count <- sum(family=="poisson")
  if(!is.null(offset)&&(poisson_count>0)) {
    if(is.vector(offset)) {
      offset <- matrix(offset, data_size, poisson_count)
    } else {
      offset <- as.matrix(offset)
    }
  }

  if (!is.null(offset))
    loffset <- log(offset)

  # check part counts
  if(length(formula)[2] < 1+additional)
    if(additional) {
      stop("Right hand side part of formula with additional variables must have at least TWO parts!")
    } else {
      stop("Right hand side part of formula must have at least ONE part!")
    }

  # theme count
  theme_R <- length(formula)[2] - additional

  # check H (number of components to keep per theme)
  H <- as.integer(H)
  if(length(H) != theme_R)
    stop(sprintf("H must be a vector of R integers! (R=number of themes=%i for this call)", theme_R))
  if(any(H<0))
    stop("H must be positive integers")
  # browser()
  if(sum(H)==0 && J==0 && !additional)
    stop("H or J must contain at least one non null value")

  # extract right-hand side (X)
  # X is a list of R formulae of the form ~...
  # As a convenience we also keep value of current 'r' and 'Hr' as attributes
  theme_X <- lapply(1:theme_R, function(r) structure(stats::terms(formula, lhs=0, rhs=r),
                                                     oldr=r, label=sprintf("T%s",r)))

  # name themes
  theme_labels <- sprintf("T%s", 1:theme_R)

  # removes themes when no components are requested (Hr=0)
  theme_X[which(H==0)] <- NULL
  theme_labels <- theme_labels[which(H>0)]
  H <- H[which(H>0)]
  theme_R <- length(H)
  if(theme_R > 0){
    for(r in 1:theme_R) {
      attr(theme_X[[r]],"r") <- r
      attr(theme_X[[r]],"Hr") <- H[r]
    }
  }


  # extract additional variables (A)
  # A is a formula of the form ~...
  if(additional) {
    terms_A <- stats::terms(formula, lhs=0, rhs=length(formula)[[2]])
  } else {
    terms_A <- NULL
  }

  # extract vars
  data_vars <- names(data)
  Theme_vars <- lapply(theme_X, all.vars)
  X_vars <- unique(unlist(Theme_vars))
  A_vars <- all.vars(terms_A)
  YXA_vars <- unique(c(Y_vars, X_vars, A_vars))

  # check if all variables can be found within data
  missing_vars <- YXA_vars[!YXA_vars %in% data_vars]
  if(length(missing_vars))
    stop("Some variable(s) where not found in data! '", paste(missing_vars, collapse="', '"),"'")

  # check that Y and X+A variable do not overlap
  error_vars <- Y_vars[Y_vars %in% c(X_vars, A_vars)]
  if(length(error_vars))
    stop("LHS and RHS variables must be different! '", paste(error_vars, collapse="', '"),"'")

  # check that A variables do not overlap with X
  error_vars <- A_vars[A_vars %in% X_vars]
  if(length(error_vars))
    stop("Additional variables must be different of theme variables! '", paste(error_vars, collapse="', '"),"'")

  # build data
  data <- data[,YXA_vars]

  # check if the number of factors is lower than the number of responses
  J.max <- choix.num.param(length(Y_vars))
  if(J < 0 | J > J.max)
    stop('J must be between 0 and ', J.max, "")

  #***************#
  #initialisation #----
  #***************#
  Y <- as.matrix(data[,Y_vars])
  X <- as.matrix(data[,X_vars])
  n <- nrow(Y)
  X <- apply(X, 2, FUN =  function(x) return(wtScale(x=x, w=1/n)) )
  Theme <- lapply(Theme_vars, function(c) as.matrix(X[,c]))
  if(!is.null(terms_A)) A <- as.matrix(data[,A_vars]) else A <- NULL
  K <- length(Y_vars)
  p <- ncol(X)
  if(!is.null(terms_A)) r1 <- ncol(A) else r1 <- 0
  muinf <- 1e-5

  #initialisation des u et f

  if(theme_R > 0){
    U <- list()
    F <- cbind()
    names_comp <- c()
    names_u <- c()

    for(r in 1:theme_R){
      U[[r]] <- cbind(plsr(formula = Y~Theme[[r]], ncomp = H[r])$loading.weights)
      colnames(U[[r]]) <- paste(attr(theme_X[[r]], "label"), "_", "U", 1:H[r], sep="")
      rownames(U[[r]]) <- Theme_vars[[r]]
      F <- cbind(F, Theme[[r]]%*%U[[r]])
      names_comp <- c(names_comp, paste(attr(theme_X[[r]], "label"), "_", "SC", 1:H[r], sep=""))
    }
    colnames(F) <- names_comp
    f_old <- F
  } else {
    F <- NULL
    names_comp <- NULL
  }

  ### Initialization, Z working variables
  mu0 <- apply(Y, 2, function(x) mean(x))
  mu0 <- matrix(mu0, n, K, byrow = TRUE)
  Z <- Y

  if ("bernoulli" %in% family) {
    tmu0 <- mu0[, family == "bernoulli"]
    Z[, family == "bernoulli"] <- log(tmu0/(1 - tmu0)) +
      (Y[, family == "bernoulli"] - tmu0)/(tmu0 * (1 - tmu0))
  }
  if ("binomial" %in% family) {
    tmu0 <- mu0[, family == "binomial"]
    Z[, family == "binomial"] <- log(tmu0/(1 - tmu0)) +
      (Y[, family == "binomial"] - tmu0)/(tmu0 * (1 - tmu0))
  }
  if ("poisson" %in% family) {
    tmu0 <- mu0[, family == "poisson"]
    if (is.null(offset)) {
      Z[, family == "poisson"] <- log(tmu0) + (Y[, family == "poisson"] - tmu0)/tmu0
    } else {
      Z[, family == "poisson"] <- log(tmu0) - loffset + (Y[, family == "poisson"] - tmu0)/tmu0
    }
  }
  if ("gaussian" %in% family) {
    Z[, family == "gaussian"] <- Y[, family == "gaussian"]
  }

  # initialisation des eta

  if (is.null(A)) {
    reg <- cbind(1, F)
    sol <- solve(crossprod(reg), crossprod(reg, Z))
    eta.old <- reg%*%sol
  } else {
    reg <- cbind(1, A, F)
    sol <- solve(crossprod(reg), crossprod(reg, Z))
    eta.old <- reg%*%sol
  }

  # Update initialization of Z and initialization of W Z <- eta
  sigma2.old <- NULL
  W <- matrix(1, n, K)

  if ("bernoulli" %in% family) {
    etainf <- log(muinf/(1 - muinf))
    indinf <- 1 * (eta.old[, family == "bernoulli"] < etainf)
    eta.old[, family == "bernoulli"] <- eta.old[, family == "bernoulli"] * (1 - indinf) + etainf * indinf
    indsup <- 1 * (eta.old[, family == "bernoulli"] > -etainf)
    eta.old[, family == "bernoulli"] <- eta.old[, family == "bernoulli"] * (1 - indsup) - etainf * indsup
    mu <- exp(eta.old[, family == "bernoulli"])/(1 + exp(eta.old[, family == "bernoulli"]))
    Z[, family == "bernoulli"] <- eta.old[, family == "bernoulli"] + (Y[, family == "bernoulli"] - mu)/(mu *
                                                                                                          (1 - mu))
    W[, family == "bernoulli"] <- mu * (1 - mu)
  }
  if ("binomial" %in% family) {
    etainf <- log(muinf/(1 - muinf))
    indinf <- 1 * (eta.old[, family == "binomial"] < etainf)
    eta.old[, family == "binomial"] <- eta.old[, family == "binomial"] * (1 - indinf) + etainf * indinf
    indsup <- 1 * (eta.old[, family == "binomial"] > -etainf)
    eta.old[, family == "binomial"] <- eta.old[, family == "binomial"] * (1 - indsup) - etainf * indsup
    mu <- exp(eta.old[, family == "binomial"])/(1 + exp(eta.old[, family == "binomial"]))
    Z[, family == "binomial"] <- eta.old[, family == "binomial"] + (Y[, family == "binomial"] - mu)/(mu *
                                                                                                       (1 - mu))
    W[, family == "binomial"] <- mu * (1 - mu) * size
  }
  if ("poisson" %in% family) {
    etainf <- log(muinf)
    indinf <- 1 * (eta.old[, family == "poisson"] < etainf)
    eta.old[, family == "poisson"] <- eta.old[, family == "poisson"] * (1 - indinf) + etainf * indinf
    indsup <- 1 * (eta.old[, family == "poisson"] > -etainf)
    eta.old[, family == "poisson"] <- eta.old[, family == "poisson"] * (1 - indsup) - etainf * indsup
    if (is.null(offset)) {
      mu <- exp(eta.old[, family == "poisson"])
    } else {
      mu <- exp(eta.old[, family == "poisson"] + loffset)
    }
    Z[, family == "poisson"] <- eta.old[, family == "poisson"] + (Y[, family == "poisson"] - mu)/mu
    W[, family == "poisson"] <- mu
  }
  if ("gaussian" %in% family) {
    Z[, family == "gaussian"] <- Y[, family == "gaussian"]
    sigma2.old <- rep(0.1, sum(family == "gaussian"))
    ind.gauss <- which(family == "gaussian")
    for(k in 1:sum(family=="gaussian")) W[,ind.gauss[k]] <- 1/sigma2.old[k]
  } else {
    sigma2.new <- NULL
  }
  #apply(W, 2, function(x) x/sum(x))

  if(J > 0){
    B.old <- matrix(0.1, J, K)
    B.old[lower.tri(B.old)] <- 0
    G <- matrix(data=0.1,n,J)
    vc.old <- crossprod(B.old)
  } else {
    B.old <- 0
    vc.old <- crossprod(B.old)
    G <- 0
  }


  #***********#
  # main loop #----
  #***********#
  if("gaussian" %in% family) tol <- rep(Inf, 5)
  else tol <- rep(Inf, 4)
  a <- 1
  while(any(tol > method$epsilon) & a < method$maxiter){

    ## EM----

    if(J > 0){
      res.EM <- EM.glm.factorSCGLR(B=B.old, sol=sol, G=G, sigma2=sigma2.old, crit=crit,
                                   W=W, J=J, A=A, F=F, Z=Z, family=family)


      B.new <- res.EM$B
      sigma2.new <- res.EM$sigma2
      sol <- res.EM$sol
      G <- res.EM$G
      vc.new <- crossprod(B.new)
    } else {
      B.new <- 0
      G <- 0
      vc.new <- crossprod(B.new)
    }


    ## FSA----
    if(J > 0){
      if (is.null(A)) {
        reg <- cbind(1, F)
        eta.new <- reg%*%sol + G%*%B.new
      } else {
        reg <- cbind(1, A, F)
        eta.new <- reg%*%sol + G%*%B.new
      }
    } else {
      if (is.null(A)) {
        reg <- cbind(1, F)
        # browser()
        for (j in seq(K)) {
          sol[,j] <- solve(crossprod(reg, W[, j] * reg),
                           crossprod(reg, W[, j] * Z[,j]))
        }
        eta.new <- reg%*%sol
      } else {
        reg <- cbind(1, A, F)
        for (j in seq(K)) {
          sol[,j] <- solve(crossprod(reg, W[, j] * reg),
                           crossprod(reg, W[, j] * Z[,j]))
        }
        eta.new <- reg%*%sol
      }
      if("gaussian" %in% family){
        sigma2.new <- colMeans((Z[, family == "gaussian"]-eta.new[, family == "gaussian"])^2)
      }
    }

    # browser()

    if ("bernoulli" %in% family) {
      etainf <- log(muinf/(1 - muinf))
      indinf <- 1 * (eta.new[, family == "bernoulli"] < etainf)
      eta.new[, family == "bernoulli"] <- eta.new[, family == "bernoulli"] * (1 - indinf) + etainf * indinf
      indsup <- 1 * (eta.new[, family == "bernoulli"] > -etainf)
      eta.new[, family == "bernoulli"] <- eta.new[, family == "bernoulli"] * (1 - indsup) - etainf * indsup
      mu <- exp(eta.new[, family == "bernoulli"])/(1 + exp(eta.new[, family == "bernoulli"]))
      Z[, family == "bernoulli"] <- eta.new[, family == "bernoulli"] + (Y[, family == "bernoulli"] - mu)/(mu *
                                                                                                            (1 - mu))
      W[, family == "bernoulli"] <- mu * (1 - mu)
    }
    if ("binomial" %in% family) {
      etainf <- log(muinf/(1 - muinf))
      indinf <- 1 * (eta.new[, family == "binomial"] < etainf)
      eta.new[, family == "binomial"] <- eta.new[, family == "binomial"] * (1 - indinf) + etainf * indinf
      indsup <- 1 * (eta.new[, family == "binomial"] > -etainf)
      eta.new[, family == "binomial"] <- eta.new[, family == "binomial"] * (1 - indsup) - etainf * indsup
      mu <- exp(eta.new[, family == "binomial"])/(1 + exp(eta.new[, family == "binomial"]))
      Z[, family == "binomial"] <- eta.new[, family == "binomial"] + (Y[, family == "binomial"] - mu)/(mu *
                                                                                                         (1 - mu))
      W[, family == "binomial"] <- mu * (1 - mu) * size
    }
    if ("poisson" %in% family) {
      etainf <- log(muinf)
      indinf <- 1 * (eta.new[, family == "poisson"] < etainf)
      eta.new[, family == "poisson"] <- eta.new[, family == "poisson"] * (1 - indinf) + etainf * indinf
      indsup <- 1 * (eta.new[, family == "poisson"] > -etainf)
      eta.new[, family == "poisson"] <- eta.new[, family == "poisson"] * (1 - indsup) - etainf * indsup
      if (is.null(offset)) {
        mu <- exp(eta.new[, family == "poisson"])
      } else {
        mu <- exp(eta.new[, family == "poisson"] + loffset)
      }
      Z[, family == "poisson"] <- eta.new[, family == "poisson"] + (Y[, family == "poisson"] - mu)/mu
      W[, family == "poisson"] <- mu
    }
    if ("gaussian" %in% family) {
      Z[, family == "gaussian"] <- Y[, family == "gaussian"]
      for(k in 1:sum(family=="gaussian")) W[,ind.gauss[k]] <- 1/sigma2.new[k]
    }


    ## PING----
    if(theme_R > 0){
      for(r in 1:theme_R){
        names_r <- paste(attr(theme_X[[r]], "label"), "_", "SC", 1:H[r], sep="")

        if(theme_R > 1 | additional) Ar <- cbind(A, F[,!names_comp%in%names_r])
        else Ar <- NULL
        for(h in 1:H[r]){
          if(h==1){
            U[[r]][,h] <- ping(Z=Z,X=Theme[[r]],AX=Ar,W=W,F=NULL,
                               u=U[[r]][,h],method=method)
          } else {
            u <- U[[r]][,h]
            C <- crossprod(Theme[[r]], F[,names_r[1:(h-1)]])
            u <- u - C %*% solve(crossprod(C), crossprod(C,u))
            u <- c(u/sqrt(sum(u^2)))
            U[[r]][,h] <- ping(Z=Z,X=Theme[[r]],AX=Ar,W=W,F=F[,names_r[1:(h-1)]],
                               u=u,method=method)
          }
          F[,names_r[h]] <- cbind(Theme[[r]]%*%U[[r]][,h])
        }
      }
      f_old <- apply(f_old,2,function(x) x/(sqrt(c(crossprod(x)/n))))
      f_new <- apply(F,2,function(x) x/(sqrt(c(crossprod(x)/n))))
      tol[1] <- abs(sum(1 - diag(crossprod(f_old, f_new)/n)^2))
      f_old <- F
    } else {
      tol[1] <- 0
    }

    tol[2] <- mean((B.old-B.new)^2)
    B.old <- B.new
    tol[3] <- mean((eta.old-eta.new)^2)
    eta.old <- eta.new
    tol[4] <- mean((vc.old-vc.new)^2)
    vc.old <- vc.new
    if(!is.null(sigma2.old)){
      tol[5] <- mean((sigma2.old-sigma2.new)^2)
      sigma2.old <- sigma2.new
    }

    a <- a+1
  }# end main loop

  #******#
  # out  #----
  #******#

  coef <- sol
  names.coef <- c("intercept", A_vars, names_comp)
  row.names(coef) <- names.coef
  if(theme_R == 0){
    Theme <- NULL
    U <- NULL
  }
  if(J > 0){
    if(J > 1){
      rotation <- rep(1, J)
      for(k in 1:J) if(B.new[k,k]<0) rotation[k] <- -1
      B.rotation <- diag(rotation) %*% B.new
      G.rotation <- res.EM$G %*% diag(rotation)
      rownames(B.rotation) <- paste("Factor", 1:J, sep="")
      colnames(G.rotation) <- paste("Factor", 1:J, sep="")
    } else {
      if(B.new[1,1]<0){
        B.rotation <- -B.new
        G.rotation <- -res.EM$G
        rownames(B.rotation) <- paste("Factor", 1:J, sep="")
        colnames(G.rotation) <- paste("Factor", 1:J, sep="")
      }
    }

  } else {
    G.rotation <- NULL
    B.rotation <- NULL
  }

  out <- list(U=U,
              comp=F,
              B=B.rotation,
              sigma2=sigma2.new,
              coef=coef,
              eta=eta.new,
              G=G.rotation,
              offset=offset,
              Z=Z,
              Theme=Theme,
              Y=Y,
              A=A,
              W=W,
              family=family)

  class(out) <- "FactorsSCGLR"

  return(out)
}





EM.glm.factorSCGLR <- function(B, sigma2, sol, G, J, A, W, F, Z, crit, family){
  K <- ncol(Z)
  N <- nrow(Z)
  G.old <- G
  G.new <- matrix(data=0, nrow = nrow(G), ncol=ncol(G))
  B.old <- B
  B.new <- matrix(data=0, ncol = ncol(B), nrow = nrow(B))
  if(!is.null(sigma2)){
    sigma2.old <- sigma2
    sigma2.new <- rep(0, length(sigma2))
  } else {
    sigma2.new <- NULL
  }
  sol.old <- sol
  sol.new <- matrix(data=0, ncol = ncol(sol), nrow = nrow(sol))
  vc.old <- crossprod(B.old)
  if(!is.null(sigma2)) tol <- rep(Inf, 5)
  else tol <- rep(Inf, 4)
  if (is.null(A)) reg <- cbind(1, F)
  else reg <- cbind(1, A, F)

  i <- 0
  # browser()
  while( any(tol > crit$tol) & i < crit$maxit){

    #*******#
    #E step # ----
    #*******#
    R <- list()

    BB <- crossprod(B.old)
    for(n in 1:N){
      alpha <- B.old%*%spdinv(BB+Diag.matrix(K,1/W[n,]))
      G.new[n,] <- alpha%*%(Z[n,]-crossprod(sol.old, reg[n,]))
      G.new[n, G.new[n,] > 1e5] <- 1e5
      G.new[n, G.new[n,] < -1e5] <- -1e5
      R[[n]] <- Diag.matrix(J,1)-tcrossprod(alpha, B.old)+tcrossprod(G.new[n,])
    }

    #*******#
    #M step # ----
    #*******#

    # actualisation des parametres de regression
    Z_bar <- Z - G.new %*% B.old
    for (j in seq(K)) {
      sol.new[,j] <- solve(crossprod(reg, W[, j] * reg),
                           crossprod(reg, W[, j] * Z_bar[,j]))
    }


    # actualisation des sigma2
    Z_tild <- Z-reg%*%sol.new

    if(!is.null(sigma2)){
      R_tild <- Reduce(`+`, R)
      ind.gauss <- which(family == "gaussian")
      for(k in 1:sum(family=="gaussian")){
        sigma2.new[k] <- (crossprod(Z_tild[,ind.gauss[k]])+
                            crossprod(B.old[,ind.gauss[k]], R_tild%*%B.old[,ind.gauss[k]])-
                            2*crossprod(G.new%*%B.old[,ind.gauss[k]], Z_tild[,ind.gauss[k]]))/N
      }
    }

    #actualisation des b
    for(k in 1:K){
      R1 <- matrix(data=0, nrow=J, ncol=J)
      for(n in 1:N) R1 <- R1 + W[n,k] * R[[n]]

      if(k < J+1){
        B.new[1:k,k] <- solve(R1[1:k,1:k], crossprod(G.new[,1:k], W[, k] * Z_tild[,k]))
      } else {
        B.new[,k] <- solve(R1, crossprod(G.new, W[, k] * Z_tild[,k]))
      }
    }
    B.new[B.new > 1e5] <- 1e5
    B.new[B.new < -1e5] <- -1e5

    vc.new <- crossprod(B.new)

    i <- i+1

    tol[1] <- mean((B.old-B.new)^2)
    tol[2] <- mean((sol.old-sol.new)^2)
    tol[3] <- mean((vc.old-vc.new)^2)
    tol[4] <- mean((G.old-G.new)^2)
    if(!is.null(sigma2)) tol[5] <- mean((sigma2.old-sigma2.new)^2)


    B.old <- B.new
    if(!is.null(sigma2)) sigma2.old <- sigma2.new
    vc.old <- vc.new
    sol.old <- sol.new
    G.old <- G.new

  }# end loop EM

  return(list(B=B.new, sigma2=sigma2.new, sol=sol.new, G=G.new))

}
