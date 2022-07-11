if(getRversion()>="2.15.1") {
  # remove warnings due to ggplot2 syntax strangeness
  utils::globalVariables(c("x","y","label","angle","hjust"))
}

#' @title FactorSCGLR generic plot
#' @param x an object from FactorSCGLR
#' @param thresold correlations with the two components of the plane lower than this threshold will be ignored
#' @param theme a integer indicating the theme
#' @param plan a size-2 vector indicating which components are plotted
#' @param lin.pred a logical. Should we draw linear predictor
#' @return an object of class ggplot2
#' @export
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
#'res <- FactorSCGLR(formula=form, data=data, H=H, J=J,
#'                    family=fam, method=met, offset = data$surface)
#'
#'# show the correlation plot
#'plot_comp(x=res, thresold=0.5, theme=1, plan=c(1,2), lin.pred=TRUE)

plot_comp <- function(x, thresold=0, plan=c(1,2), theme=1, lin.pred=FALSE){
  if (class(x) != "FactorSCGLR")
    stop("This plot function need an FactorSCGLR result")

  res <- x
  labels.offset <- 0.01

  n_comp <- colnames(res$comp)
  number <- paste("T", theme, sep = "")
  comps <- n_comp[grep(number, n_comp)]

  n_var <- colnames(res$Theme)
  var <- n_var[grep(number, n_var)]

  inertia <- cor(res$Theme[, var], res$comp[, comps])
  inertia <- inertia^2
  inertia <- colMeans(inertia)
  # browser()
  p <- qplot()+
    coord_fixed()+ theme(axis.title = element_text(size = 30),
                         plot.title = element_text(size = 30)) +
    labs(title = paste("Component plane (", plan[1], ",", plan[2], ")",
                       " for Theme ", theme, sep = "")) +
    # thicker x unit arrow
    xlab(paste("SC", plan[1], " (", round(100*inertia[plan[1]],2), "%", ")", sep = "")) +
    geom_hline(yintercept=0)+
    geom_segment(aes(x=-1.1,xend=1.1,y=0,yend=0),size=1,arrow=arrow(length=unit(0.02,"npc")))+
    # thicker y unit arrow
    ylab(paste("SC", plan[2], " (", round(100*inertia[plan[2]],2), "%", ")", sep = "")) + #theme(axis.title.y = element_text(size = 20)) +
    geom_vline(xintercept=0)+
    geom_segment(aes(y=-1.1,yend=1.1,x=0,xend=0),size=1,arrow=arrow(length=unit(0.02,"npc")))
  p <- p + annotation_custom(circleGrob(r=0.5,gp=gpar(fill=NA)),-1,1,-1,1)
  p <- p + annotation_custom(circleGrob(r=0.5,gp=gpar(lty=2, fill=NA)),
                             -thresold,thresold,-thresold,thresold)

  co1 <- as.data.frame(cor(res$Theme[, var],
                           res$comp[,comps][,c(plan[1], plan[2])]))
  names(co1) <- c("x", "y")
  co1$norm <- sqrt(co1$x^2+co1$y^2)
  co1$label <- names(as.data.frame(res$Theme[, var]))
  co1$arrows.color <- 'black'
  co1$labels.color <- 'black'
  co1$labels.size <- 6
  co1$angle <- atan2(co1$y,co1$x)*180/pi
  co1$hjust <- ifelse(abs(co1$angle)>90,1,0)
  co1$angle <- ifelse(abs(co1$angle)>90,co1$angle+180,co1$angle)

  co1 <- co1[co1$norm>thresold,]

  p <- p + geom_segment(
    aes(x=0,y=0,xend=x,yend=y),
    data=co1,
    color=co1$arrows.color,
    arrow=arrow(length=unit(0.02,"npc"))
  )

  p <- p + geom_text(
    aes(x=x*(1+labels.offset/norm),y=y*(1+labels.offset/norm),label=label,angle=angle,hjust=hjust),
    data=co1,
    color=co1$labels.color,
    size=co1$labels.size
  )
  # browser()
  if(lin.pred==TRUE){
    co2 <- as.data.frame(cor(cbind(1, res$A, res$comp) %*% res$coef,
                             res$comp[,comps][,c(plan[1], plan[2])]   ))
    # browser()
    names(co2) <- c("x", "y")
    co2$norm <- sqrt(co2$x^2+co2$y^2)
    co2$label <- names(as.data.frame(res$Y))
    co2$arrows.color <- 'red'
    co2$labels.color <- 'red'
    co2$labels.size <- 6
    co2$angle <- atan2(co2$y,co2$x)*180/pi
    co2$hjust <- ifelse(abs(co2$angle)>90,1,0)
    co2$angle <- ifelse(abs(co2$angle)>90,co2$angle+180,co2$angle)
    co2$alpha <- 0.3

    co2 <- co2[co2$norm>thresold,]

    p <- p + geom_segment(
      aes(x=0,y=0,xend=x,yend=y),
      data=co2,
      color=co2$arrows.color,
      alpha=co2$alpha,
      arrow=arrow(length=unit(0.02,"npc"))
    )

    p <- p + geom_text(
      aes(x=x*(1+labels.offset/norm),y=y*(1+labels.offset/norm),label=label,angle=angle,hjust=hjust),
      data=co2,
      color=co2$labels.color,
      alpha=co2$alpha,
      size=co2$labels.size
    )

  }

  return(list(p=p))
}
