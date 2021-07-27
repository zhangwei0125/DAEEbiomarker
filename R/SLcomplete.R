#' Efficient estimation for biomarker studies embedded in randomized clinical trials
#' in the complete-data setting.
#'  
#' 'SLcomplete()' returns efficient estimation for biomarker studies in the complete-
#' data setting in which the biomarker is measured for all subjects in the trial, by 
#' using machine learning method.
#'   
#' This function is used to calculate estimates for the relationship
#' between a clinical outcome of interest and the biomarker, treatment indicator (can
#' be omitted in the case of a prognostic biomarker) in the complete-data setting, by 
#' using the proposed efficient estimation approach. The efficient estimation incorporates
#' the independence between treatment and all relevant baseline covariates. 
#'   
#' @param Y a vector for outcome
#' @param Trt a vector for treatment indicator
#' @param Z a vector for biomarker
#' @param W a matrix for baseline covariates (not inlcuing the biomarker)
#' @param Pi probability of receiving the treatment (T=1) 
#' @param SL.lib a list of functions for candidate prediction algorithms
#' @param SL.family allows gaussian or binomial to describe the error distribution
#' @param SL.method a character value indicating the method used to estimate the coefficients
#' for the super learner and the model to combine the individual algorithms
#' @param SL.cv  Number of splits for the cross-validation step in super learner
#' @param Split.cv Number of splits for the calculation of cross-validated variance estimate
#' @param ols.max upper bound for the coefficients of individual algorithms in the super learner
#' 
#' @return A list with the first argument being the estimate of the association (beta) between the
#'  outcome and biomarker, treatment indicator for the individual algorithms, as well as a super learner
#'  that combines the individual algorithms, the second argument being the corresponding cross-validated
#'  variance estimate.
#'  
#'  @import SuperLearner rootSolve
#'  
#'  @export
#'  
SLcomplete <- function(Y, Trt, Z, W, Pi, SL.lib, SL.family, SL.method, SL.cv, Split.cv, ols.max){
  
  n <- length(Y)
  Q <- 4
  if(is.vector(W))  W <- matrix(W, ncol=1)
  datMat <- cbind(Y, Trt, Z)
  K <- length(SL.lib)
  ############################ Point estimate of beta #######################
  ################# Estimate a() based on the initial estimate beta.tilde 
  a.hat <- array(NA, dim=c(n, Q, K+1))
  a.hat[,,1:K] <- augFun.Est(Y=Y, Trt=Trt, Z=Z, W=W, newZ=Z, newW=W, Pi=Pi,
                             SL.lib=SL.lib, SL.family=SL.family)$a.fun
  
  ### calculate the coefficients for super learner estimator
  SL.coef <- calWeight(Y=Y, Trt=Trt, Z=Z, W=W, Pi=Pi, SL.lib=SL.lib, SL.family=SL.family,
                       SL.method=SL.method, SL.cv=SL.cv, ols.max=ols.max)
  for(q in 1:Q)
  {
    a.hat[,q,K+1] <- a.hat[,q,1:K] %*% matrix(SL.coef[q,], ncol=1)
  }
  
  ################# Estimate beta based on augmentated estimation equation psi.plus
  bet.plus <- matrix(NA, nrow=Q, ncol=K+1)
  var.plus <- array(NA, dim=c(Q, Q, K+1))
  for(k in 1:(K+1))
  {
    #### augmentated estimating equation
    psi.plus.Equ <- function(bet)
    {
      temp.psi <- t(apply(datMat, MARGIN=1, FUN=psi.fun, bet=bet))
      psi.plus <- temp.psi-matrix(Trt-Pi, nrow=n, ncol=Q, byrow=F)*a.hat[,,k]
      return(colSums(psi.plus))
    }
    bet.plus[,k] <- multiroot(f=psi.plus.Equ, start=rep(0, Q), maxiter=500)$root
    ##### the varaince based on the  estimate of beta and a()
    var.plus[,,k] <- calVar(Y=Y, Trt=Trt, Z=Z, Pi=Pi, bet=bet.plus[,k], a.fun=a.hat[,,k])
  }
  
  
  ##### the varaince based on the cross-validated estimate of beta and a()
  a.hat.SS <- array(NA, dim=c(n, Q, K+1))
  fold.SS <- sample(1:Split.cv, n, replace=TRUE)
  for(h in 1:Split.cv)
  {
    SS.train <- which(fold.SS!=h)
    SS.valid <- which(fold.SS==h)
    Y.SS <- Y[SS.train]; Trt.SS <- Trt[SS.train]; Z.SS <- Z[SS.train]; W.SS <- W[SS.train,]
    newZ.SS <- Z[SS.valid]; newW.SS <- W[SS.valid,]
    temp.SS <- augFun.Est(Y=Y.SS, Trt=Trt.SS, Z=Z.SS, W=W.SS, newZ=newZ.SS,
                          newW=newW.SS, Pi=Pi, SL.lib=SL.lib, SL.family=SL.family)
    a.hat.SS[SS.valid,,1:K] <- temp.SS$a.fun
    SL.coef.SS <- calWeight(Y=Y.SS, Trt=Trt.SS, Z=Z.SS, W=W.SS, Pi=Pi, SL.lib=SL.lib, SL.family=SL.family, 
                            SL.method=SL.method, SL.cv=SL.cv, ols.max=ols.max)
    for(q in 1:Q)
    {
      a.hat.SS[SS.valid,q,K+1] <- a.hat.SS[SS.valid,q,1:K] %*% matrix(SL.coef.SS[q,], ncol=1)
    }
  }
  ### update the estimate of beta
  bet.plus.SS <- matrix(NA, nrow=Q, ncol=K+1)
  var.plus.SS <- array(NA, dim=c(Q, Q, K+1))
  for(k in 1:(K+1))
  {
    #### augmentated estimating equation
    psi.plus.Equ.SS <- function(bet)
    {
      temp.psi <- t(apply(datMat, MARGIN=1, FUN=psi.fun, bet=bet))
      psi.plus.SS <- temp.psi-matrix(Trt-Pi, nrow=n, ncol=Q, byrow=F)*a.hat.SS[,,k]
      return(colSums(psi.plus.SS))
    }
    bet.plus.SS[,k] <- multiroot(f=psi.plus.Equ.SS, start=rep(0, Q), maxiter=500)$root
    var.plus.SS[,,k] <- calVar(Y=Y, Trt=Trt, Z=Z, Pi=Pi, bet=bet.plus.SS[,k], a.fun=a.hat.SS[,,k])
  }
  
  ################################################################
  ### return the results (point estimate and variance estimate of beta)
  return(list(pEst=bet.plus, varEst.cv=var.plus.SS))
}

