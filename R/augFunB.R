#' Estimation of the second agumented term in the two-phase design by indiivdual
#'  machine learning prediction algorithms
#' 
#' 'augFunB()' returns the estimate of second agumented function b(Y, T, W) in the two-phase
#'  design by indiivdual machine learning prediction algorithms
#' 
#' This function is used to estimate the second agumented term of the estimation function 
#' for the two-phase design. This agumented function is a function of (outcome,treatment, 
#' baseline covariate).
#' 
#' @param Y a vector for outcome
#' @param Trt a vector for treatment indicator
#' @param Z a vector for biomarker
#' @param R an indicator for whether the biomarker is observed 
#' @param W a matrix for baseline covariates (not inlcuing the biomarker)
#' @param Pi probability of receiving the treatment (T=1) 
#' @param probR conditional probability for the event that biomarker is observed
#' @param newY a new vector for outcome at which the agumented function b() is evaluated 
#' @param newTrt a new vector for treatment indicator at which the agumented function b() is evaluated
#' @param newW a new matrix for baseline covariates at which the agumented function b() is evaluated 
#' @param SL.lib a list of functions for candidate prediction algorithms
#' @param SL.family allows gaussian or binomial to describe the error distribution
#' @param a.fun the first augmented term which is a function of (biomarker, baseline covairate)
#' 
#' @return a list with the first argument being the initial estimate of beta without augmentation and
#' the second argument being the estimate of the augmented function b() at (newY, newTrt, newW)
#' 
#' @export
#' 
augFunB <- function(Y, Trt, Z, R, W, Pi, probR, newY, newTrt, newW, SL.lib, SL.family, a.fun){
  
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
  n.run <- length(Y)
  n.test <- length(newY)
  H <- as.data.frame(rbind(cbind(Y, Trt, W), cbind(newY, newTrt, newW)))
  colnames(H) <- paste("v", 1:ncol(H), sep="")
  H.run <- H[1:n.run,]
  H.test <- H[(1+n.run):(n.run+n.test),]
  b.lib <- array(NA, dim=c(n.test, Q, K))
  for(q in 1:Q)
  {
    ### use the b.old based on the super learner algorithm to iterate a
    rsp.b <- psi.Ini[,q]- (Trt-Pi)*a.fun[,q]
    wt.b <-  R*(1-probR)/(probR^2)
    ###
    for(k in 1:K)
    {
      b.lib[,q,k] <- SL.lib[[k]](Y=rsp.b, X=H.run, newX=H.test, family=SL.family, obsWeights=wt.b)$pred
    }
  }
  ################################################################
  return(list(bet.Ini=bet.Ini, lib.predict=b.lib))
}