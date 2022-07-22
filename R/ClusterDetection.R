#' @title Clustering steps to detect responses with mutual dependencies
#' @param mat a matrice of loadings of size J (number of factors) times K (number of responses)
#' @return \item{cluster}{a vector of integers indicating the cluster to which each response is allocated}
#' @return \item{mds}{a matrix whose rows give the coordinates of the points chosen to represent the dissimilarities}
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
#'# Detect the clusters
#'cluster <- ClusterDetection(res$B)



ClusterDetection <- function(mat){

  # compute the conditional variance covariance matrix
  covariance <- crossprod(mat)

  # compute the conditional correlation matrix
  correlation <- matrix(0, ncol(covariance), ncol(covariance))
  for(i in 1:ncol(covariance)){
    for(j in 1:ncol(covariance)){
      correlation[i,j] <- covariance[i,j]/(sqrt(covariance[i,i])*sqrt(covariance[j,j]))
    }
  }

  # compute the dissimilarity matrix
  dist <- 2*abs(1-correlation*correlation)

  # perform multidimensional scaling
  flag <- TRUE
  nb_axes <- 1
  while(flag){
    nb_axes <- nb_axes+1
    res.CMD <- cmdscale(dist,eig=TRUE, k=nb_axes)
    if(res.CMD$GOF[2] > 0.99) flag <- FALSE
  }

  # perform K-means with CAH as initialization
  k.opt <- which.max(fviz_nbclust(res.CMD$points, FUNcluster = hkmeans,
                                  method = "silhouette",
                                  k.max=ncol(mat)-1)$data$y)
  res.km <- hkmeans(x=res.CMD$points, k = k.opt)

  #out
  mds <- res.CMD$points
  rownames(mds) <- colnames(mat)
  colnames(mds) <- paste("Dim", 1:ncol(mds), sep="")

  cluster <- as.data.frame(res.km$cluster)
  cluster <- t(cluster)
  colnames(cluster) <- colnames(mat)
  rownames(cluster) <- "cluster"

  return(list(cluster = cluster,
              mds = mds))
}

