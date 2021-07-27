#'  Machine learning estimation of the second agumented term in the two-phase design
#'  using a super learner procedure
#' 
#' 'SuperLearnerB()' returns the estimate of second agumented function b() in the two-phase
#'  design using a super learner procedure
#' 
#' This function is used to estimate the second agumented term of the estimation function 
#' for the two-phase design using a super learner procedure. This agumented function is
#'  a function of (outcome, treatment, baseline covariate).
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
#' @param SL.method a character value indicating the method used to estimate the coefficients
#' for the super learner and the model to combine the individual algorithms
#' @param a.fun the first augmented term which is a function of (biomarker, baseline covairate)
#' @param ols.max upper bound for the coefficients of individual algorithms in the super learner
#' @param SL.cv  Number of splits for the cross-validation step in super learner
#' 
#' @return a list with the first argument being the estimate of augmented function b() by the individual
#' prediction algorithms and the second argument being the estimate by the super learner which combines
#' the individual algorithms.
#' 
#'  @export
#'   
SuperLearnerB <- function(Y, Trt, Z, R, W, Pi, probR, newY, newTrt, newW, SL.lib, SL.family, SL.method, a.fun, ols.max, SL.cv){
  
  n <- length(Y)
  Q <- 4
  if(is.vector(W))  W <- matrix(W, ncol=1)
  datMat <- cbind(Y, Trt, Z)
  K <- length(SL.lib)
  ######## Estimators for the algortihms in the list of library
  b.hat <- augFunB(Y=Y, Trt=Trt, Z=Z, R=R, W=W, Pi=Pi, probR=probR, newY=newY, newTrt=newTrt, newW=newW,
                   SL.lib=SL.lib, SL.family=SL.family, a.fun=a.fun)$lib.predict
  
  ######## calculate the coefficients for super learner estimator
  fold <- sample(1:SL.cv, n, replace=TRUE)
  b.hat.cv <- array(NA, dim=c(n, Q, K))
  rsp.b.cv <- matrix(NA, nrow=n, ncol=Q)
  for(v in 1:SL.cv)
  {
    train <- which(fold!=v)
    valid <- which(fold==v)
    n.valid <- length(valid)
    temp.cv <- augFunB(Y=Y[train], Trt=Trt[train], Z=Z[train], R=R[train], W=W[train,], Pi=Pi, probR=probR[train], newY=Y[valid],
                       newTrt=Trt[valid], newW=W[valid,], SL.lib=SL.lib, SL.family=SL.family, a.fun=a.fun[train,])
    b.hat.cv[valid,,] <- temp.cv$lib.predict
    bet.Ini.train <- temp.cv$bet.Ini
    datMat.valid <- cbind(Y=Y[valid], Trt=Trt[valid], Z=Z[valid])
    temp.psi.valid <- t(apply(datMat.valid, MARGIN=1, FUN=psi.fun, bet=bet.Ini.train))
    temp1.valid <- temp.psi.valid-matrix(Trt[valid]-Pi, nrow=n.valid, ncol=Q, byrow=F)*a.fun[valid,]
    temp2.valid <- matrix(sqrt(R[valid]*(1-probR[valid]))/probR[valid], nrow=n.valid, ncol=Q, byrow=F)
    rsp.b.cv[valid,] <- temp1.valid*temp2.valid
  }
  
  ### 
  res.coef <- matrix(NA, nrow=Q, ncol=K)
  for(q in 1:Q)
  {
    reg.b <- matrix(sqrt(R*(1-probR))/probR, nrow=n, ncol=K)*b.hat.cv[,q,]
    if(SL.method=="ols"){
      run.coef <- coef(lm(rsp.b.cv[,q]~0+reg.b))
      run.coef[is.na(run.coef)] <- 0
      res.coef[q,] <- pmax(pmin(run.coef, ols.max), -ols.max)
    }
    if(SL.method=="nnls"){
      run.coef <- nnls(reg.b, rsp.b.cv[,q])$x
      run.coef[is.na(run.coef)] <- 0
      if(sum(run.coef)>0)  
      {
        res.coef[q,] <- run.coef/sum(run.coef)
      }else{
        warning("All algorithms have zero weight", call.=FALSE)
        res.coef[q,] <-  run.coef
      }
    }
  }
  
  #### the super learner predictor
  b.SL <- matrix(NA, nrow=dim(b.hat)[1], ncol=Q)
  for(q in 1:Q)
  {
    b.SL[,q] <- b.hat[,q,] %*% matrix(res.coef[q,], ncol=1)
  }
  ###
  return(list(lib.predict=b.hat, SL.predict=b.SL))
}

