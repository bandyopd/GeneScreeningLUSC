#---------------------------------------------------------------------------------------------
#' ECCFIC
#' @name ech2
#' @description compute the ECCFIC between two random vectors x and y
#' 
#' @param x a numeric vector/matrix
#' @param y a numeric vector/matrix
#' @param kernel a kernel to use for x, 'gaussian' or 'distance'
#' @param sigma bandwidth of the reproducing kernel, 
#'              default is heuristic median pairwise distances of x for the 'gaussian' kernel
#' @param alpha exponent on distance in (0,2] for the 'distance' kernel
#' @param bw bandwidth of the Gaussian smoothing kernel applied on y if y is continuous,
#'           bandwidths suggested by Silverman (1986) are used unless otherwise specified.
#'           
#' @return the ECCFIC between x and y
#' @example 
#' x <- rnorm(30)
#' y <- x^2 + rnorm(30)
#' ech2(x, y)

ech2 <- function(x, y, ycont=T, kernel='gaussian', sigma='default', alpha=1, bw='default'){
  
  x <- as.matrix(x)
  n <- nrow(x)
  if(n==1){return(0)}
  
  if(kernel=='gaussian'){
    if(sigma=='default'){
      sigma <- sqrt(0.5*median(dist(x)^2))
      if(sigma==0){sigma <- 0.001}
    }
  }
  
  if(identical(x,as.matrix(y))){
    if(kernel=='gaussian'){
      K <- dnorm(as.matrix(dist(x,diag=T,upper=T)), mean=0, sd=sigma)
      return((dnorm(0, 0, sigma) - mean(K))*sqrt(2*pi)*sigma)
    }else{
      if(kernel=='distance'){
        K <- 0.5*(as.matrix(dist(x, diag=T, upper=T))^alpha)
        return(mean(K))
      }
    }
  }
  
  y <- as.matrix(y)
  p <- ncol(y)
  if(bw=='default'){
    if(p==1){bw <- 1.06*sd(y)*n^(-1/5)}else{bw <- (4/n/(p+2))^(1/(p+4))*sum(diag(cov(y)))/p}
  }
  H <- diag(n) - matrix(1, n, n)/n
  if(kernel=='gaussian'){
    K <- dnorm(as.matrix(dist(x, diag=T, upper=T)), mean=0, sd=sigma)*sqrt(2*pi)*sigma
  }else{
    if(kernel=='distance'){
      normx <- matrix(apply(x*x,1,sum)^(alpha/2),n,n)
      K <- 0.5*(normx + t(normx) - as.matrix(dist(x, diag=T, upper=T))^alpha)
    }
  }
  G <- dnorm(as.matrix(dist(y, diag=T, upper=T)), mean=0, sd=bw)
  Gstar <- t(G/rowSums(G))
  return(sum(diag(K%*%H%*%Gstar%*%t(Gstar)%*%H))/n)
  
}
#---------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------
#' ECCFIC correlation
#' @name rho_ech2
#' @description compute the ECCFIC correlation between two random vectors x and y
#' 
#' @param x a numeric vector/matrix
#' @param y a numeric vector/matrix
#' @param ... other arguments that can be passed to ech2
#' 
#' @return the ECCFIC correlation between x and y
#' #' @example 
#' x <- rnorm(30)
#' y <- x^2 + rnorm(30)
#' rho_ech2(x, y)

rho_ech2 <- function(x, y, ...){
  return(ech2(x, y, ...)/ech2(x, x, ...))
}
#---------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------
#' ESIS
#' @name esis
#' @description performs sure independence screening on right-censored data using ECCFIC
#' 
#' @param x matrix of gene expression
#' @param y observed survival time
#' @param delta failure indicator
#' @param d number of genes recruited
#' @param ... other arguments that can be passed to ech2
#' 
#' @return the vector of indices selected by ESIS, along with associated ECCFIC correlations
#' 
#' @example 
#' n <- 100
#' p <- 200
#' x <- matrix(rnorm(n*p), nrow = n)
#' e <- 0.5*rnorm(n)
#' t <- exp(apply(x[,1:3],1,sum)+e) #AFT model
#' ct <- quantile(y, probs = 0.7) #30% censoring
#' y <- ifelse(t<ct,t,ct)
#' delta <- as.numeric(t<=ct)
#' esis(x, y, delta, d=38)

library(survival)
esis <- function(x, y, delta, d, ...){
  
  x <- scale(x)
  n <- nrow(x)
  p <- ncol(x)
  xord <- x[order(y),]
  resp <- data.frame(y=y, delta=delta)
  resp <- resp[order(y),]
  fit <- survfit(Surv(resp$y,resp$delta)~1)
  s <- 1 - summary(fit,times=resp$y)$surv
  w <- rep(NA, p)
  for(j in 1:p){
    ind <- complete.cases(xord[,j])
    if(sum(ind)==0){
      w[j] <- 0
      next
    }
    Fx <- rank(xord[,j], ties.method = "max")/sum(ind)
    w[j] <- rho_ech2(s[ind], Fx[ind], ...)
  }
  ow <- order(w, decreasing = T)
  ix <- ow[1:d]
  return(list(ix = ix, rho = w[ix]))
  
}
#---------------------------------------------------------------------------------------------