#' @title Estimation of the first agumented term in the two-phase design by indiivdual
#'  machine learning prediction algorithms
#'
#' @description 'augFunA()' returns the estimate of first agumented function a() in the two-phase
#'  design by indiivdual machine learning prediction algorithms
#'
#' @details This function is used to estimate the first agumented term of the estimation function
#' for the two-phase design. This agumented function is a function of (biomarker, baseline covariate).
#'
#' @param Y a vector for outcome
#' @param Trt a vector for treatment indicator
#' @param Z a vector for biomarker
#' @param R an indicator for whether the biomarker is observed
#' @param W a matrix for baseline covariates (not inlcuing the biomarker)
#' @param Pi probability of receiving the treatment (T=1)
#' @param probR conditional probability for the event that biomarker is observed
#' @param newZ a new vector for biomarker at which the agumented function a() is evaluated
#' @param newW a new matrix for baseline covariates at which the agumented function a() is evaluated
#' @param SL.lib a list of functions for candidate prediction algorithms
#' @param SL.family allows gaussian or binomial to describe the error distribution
#' @param b.fun the second augmented term which is a function of (outcome, treatment, baseline covairate)
#'
#' @return a list with the first argument being the initial estimate of beta without augmentation and
#' the second argument being the estimate of the augmented function a() at (newZ, newW)
#'
#' @export
#'
augFunA <- function(Y, Trt, Z, R, W, Pi, probR, newZ, newW, SL.lib, SL.family, b.fun){

  n <- length(Y)
  Q <- 4
  if(is.vector(W))  W <- matrix(W, ncol=1)
  datMat <- cbind(Y, Trt, Z)
  K <- length(SL.lib)
  #### unadjsted estimator
  a.Ini <- matrix(0, nrow=n, ncol=Q)
  b.Ini <- matrix(0, nrow=n, ncol=Q)
  psiStar.Ini <- function(bet)
  {
    return(colSums(psiStar(Y=Y, Trt=Trt, Z=Z, R=R, Pi=Pi, probR=probR, bet=bet, a.fun=a.Ini, b.fun=b.Ini)))
  }
  bet.Ini <- multiroot(f=psiStar.Ini, start=rep(0, Q), maxiter=500)$root
  #### estimate a using the method in SL.lib
  psi.Ini <- t(apply(datMat, MARGIN=1, FUN=psi.fun, bet=bet.Ini))
  n.run <- length(Z)
  n.test <- length(newZ)
  X <- as.data.frame(rbind(cbind(Z, W), cbind(newZ, newW)))
  colnames(X) <- paste("v", 1:ncol(X), sep="")
  X.run <- X[1:n.run,]
  X.test <- X[(1+n.run):(n.run+n.test),]
  a.lib <- array(NA, dim=c(n.test, Q, K))
  for(q in 1:Q)
  {
    ### use the b.old based on the super learner algorithm to iterate a
    rsp.a <- (psi.Ini[,q]-(1-probR)*b.fun[,q])/(Trt-Pi)
    wt.a <-  R*(Trt-Pi)^2/(probR^2)
    ###
    for(k in 1:K)
    {
      a.lib[,q,k] <- SL.lib[[k]](Y=rsp.a, X=X.run, newX=X.test, family=SL.family, obsWeights=wt.a)$pred
    }
  }
  ################################################################
  return(list(bet.Ini=bet.Ini, lib.predict=a.lib))
}
