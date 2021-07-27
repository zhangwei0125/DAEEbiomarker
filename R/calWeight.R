#' Calculation of the weights for the super learner.
#' 
#' 'calWeight()' estimates the coefficients for the super learner and the model to combine
#'  the individual algorithms in the library.
#'  
#'  This function is used to estimate the weights for the super learner and the model to combine
#'  the individual algorithms in the library.
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
#' @param ols.max upper bound for the coefficients of individual algorithms in the super learner
#' 
#' @return a vector of the coefficients for the individual algorithms in the library
#' 
#' @export
#' 
calWeight <- function(Y, Trt, Z, W, Pi, SL.lib, SL.family, SL.method, SL.cv, ols.max){
  
  n <- length(Y)
  Q <- 4
  if(is.vector(W))  W <- matrix(W, ncol=1)
  K <- length(SL.lib)
  ### calculate the super learner predictor using the cross-validation method
  fold <- sample(1:SL.cv, n, replace=TRUE)
  psi.cv <- matrix(NA, nrow=n, ncol=Q)
  a.hat.cv <- array(NA,dim=c(n, Q, K))
  for(v in 1:SL.cv)
  {
    valid <- which(fold==v)
    train <- which(fold!=v)
    temp.cv <- augFun.Est(Y=Y[train], Trt=Trt[train], Z=Z[train], W=W[train,], newZ=Z[valid],
                          newW=W[valid,], Pi=Pi, SL.lib=SL.lib, SL.family=SL.family)
    a.hat.cv[valid,,1:K] <- temp.cv$a.fun
    bet.Ini.train <- temp.cv$bet.Ini
    datMat.valid <- cbind(Y=Y[valid], Trt=Trt[valid], Z=Z[valid])
    psi.cv[valid,] <- t(apply(datMat.valid, MARGIN=1, FUN=psi.fun, bet=bet.Ini.train))
  }
  ### calculate the coefficients for super learner estimator
  res.coef <- matrix(NA, nrow=Q, ncol=K)
  for(q in 1:Q)
  {
    reg.a <- matrix(Trt-Pi, nrow=n, ncol=K)*a.hat.cv[,q,1:K]
    if(SL.method=="ols"){
      run.coef <- coef(lm(psi.cv[,q]~0+reg.a))
      run.coef[is.na(run.coef)] <- 0
      res.coef[q,] <- pmax(pmin(run.coef, ols.max), -ols.max)
    }
    if(SL.method=="nnls"){
      run.coef <- nnls(reg.a, psi.cv[,q])$x
      run.coef[is.na(run.coef)] <- 0
      if(sum(run.coef)>0)  
      {
        res.coef[q,] <- run.coef/sum(run.coef)
      }else{
        #warning("All algorithms have zero weight", call.=FALSE)
        res.coef[q,] <-  run.coef
      }
    }
  }
  ###
  return(res.coef)
}
