#' @title Machine learning estimate for the augmented function in the complete-data setting
#'
#' @description 'augFun.Est()' returns the augmented function a() at some given observations.
#'
#' @details This function is used to estimate the augmented term in the augmentation estimation
#'   function using a library of prediction algorithms.
#'
#' @param Y a vector for outcome
#' @param Trt a vector for treatment indicator
#' @param Z a vector for biomarker
#' @param W a matrix for baseline covariates (not inlcuing the biomarker)
#' @param newZ biomarker value to be evaluated
#' @param newW baseline covariate values to be evaluated
#' @param Pi probability of receiving the treatment (T=1)
#' @param SL.lib a list of functions for candidate prediction algorithms
#' @param SL.family allows gaussian or binomial to describe the error distribution
#'
#' @return A list with the first argument being the estimate of the interested
#' parameter beta based on the unaugmented estimation function and the second argument being
#' being the augmented term a() evaluated at the given observation (newZ, newW)
#'
#'  @export
#'
augFun.Est <- function(Y, Trt, Z, W, newZ, newW, Pi, SL.lib, SL.family){

  n <- length(Y)
  Q <- 4
  datMat <- cbind(Y, Trt, Z)
  psi.Equ <- function(bet)
  {
    return(rowSums(apply(datMat, MARGIN=1, FUN=psi.fun, bet=bet)))
  }
  solver <- multiroot(f=psi.Equ, start=rep(0, Q), maxiter=500)
  bet.Ini <- solver$root
  K <- length(SL.lib)
  #### estimate a using the method in SL.lib
  psi.hat <- t(apply(datMat, MARGIN=1, FUN=psi.fun, bet=bet.Ini))
  rsp <-  psi.hat/matrix(Trt-Pi, nrow=n, ncol=Q, byrow=F) # response vector in weighted prediction
  wt <- (Trt-Pi)^2 # weight vector in weighted prediction
  n.train <- length(Z)
  n.test <- length(newZ)
  X <- as.data.frame(rbind(cbind(Z, W), cbind(newZ, newW)))
  X.train <- X[1:n.train,]
  X.test <- X[(1+n.train):(n.train+n.test),]
  a.fun <- array(NA, dim=c(n.test, Q, K))
  for(q in 1:Q)
  {
    for(k in 1:K)
    {
      a.fun[,q,k] <- SL.lib[[k]](Y=rsp[,q], X=X.train, newX=X.test, family=SL.family, obsWeights=wt)$pred
    }
  }
  ################################################################
  return(list(bet.Ini=bet.Ini, a.fun=a.fun))
}
