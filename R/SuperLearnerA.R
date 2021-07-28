#' @title Machine learning estimation of the first agumented term in the two-phase design
#'  using a super learner procedure
#'
#' @description 'SuperLearnerA()' returns the estimate of first agumented function a() in the two-phase
#'  design using a super learner procedure
#'
#' @details This function is used to estimate the first agumented term of the estimation function
#' for the two-phase design using a super learner procedure. This agumented function is
#'  a function of (biomarker, baseline covariate).
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
#' @param SL.method a character value indicating the method used to estimate the coefficients
#' for the super learner and the model to combine the individual algorithms
#' @param b.fun the second augmented term which is a function of (outcome, treatment, baseline covairate)
#' @param ols.max upper bound for the coefficients of individual algorithms in the super learner
#' @param SL.cv  Number of splits for the cross-validation step in super learner
#'
#' @return a list with the first argument being the estimate of augmented function a() by the individual
#' prediction algorithms and the second argument being the estimate by the super learner which combines
#' the individual algorithms.
#'
#'  @export
#'
SuperLearnerA <- function(Y, Trt, Z, R, W, Pi, probR, newZ, newW, SL.lib, SL.family, SL.method, b.fun, ols.max, SL.cv){

  n <- length(Y)
  Q <- 4
  if(is.vector(W))  W <- matrix(W, ncol=1)
  datMat <- cbind(Y, Trt, Z)
  K <- length(SL.lib)
  ########Estimators for the algortihms in the list of library
  a.hat <- augFunA(Y=Y, Trt=Trt, Z=Z, R=R, W=W, Pi=Pi, probR=probR, newZ=newZ, newW=newW,
                   SL.lib=SL.lib, SL.family=SL.family, b.fun=b.fun)$lib.predict

  #######calculate the coefficients for super learner estimator
  fold <- sample(1:SL.cv, n, replace=TRUE)
  a.hat.cv <- array(NA,dim=c(n, Q, K))
  rsp.a.cv <- matrix(NA, nrow=n, ncol=Q)
  for(v in 1:SL.cv)
  {
    train <- which(fold!=v)
    valid <- which(fold==v)
    n.valid <- length(valid)
    temp.cv <- augFunA(Y=Y[train], Trt=Trt[train], Z=Z[train], R=R[train], W=W[train,], Pi=Pi, probR=probR[train], newZ=Z[valid],
                       newW=W[valid,], SL.lib=SL.lib, SL.family=SL.family, b.fun=b.fun[train,])
    a.hat.cv[valid,,] <- temp.cv$lib.predict
    bet.Ini.train <- temp.cv$bet.Ini
    datMat.valid <- cbind(Y=Y[valid], Trt=Trt[valid], Z=Z[valid])
    temp.psi.valid <- t(apply(datMat.valid, MARGIN=1, FUN=psi.fun, bet=bet.Ini.train))
    temp1.valid <- temp.psi.valid-matrix(1-probR[valid], nrow=n.valid, ncol=Q, byrow=F)*b.fun[valid,]
    temp2.valid <- matrix(sqrt(R[valid])/probR[valid], nrow=n.valid, ncol=Q, byrow=F)
    rsp.a.cv[valid,] <- temp1.valid*temp2.valid
  }
  ###
  res.coef <- matrix(NA, nrow=Q, ncol=K)
  for(q in 1:Q)
  {
    reg.a <- matrix(sqrt(R)*(Trt-Pi)/probR, nrow=n, ncol=K)*a.hat.cv[,q,]
    if(SL.method=="ols"){
      run.coef <- coef(lm(rsp.a.cv[,q]~0+reg.a))
      run.coef[is.na(run.coef)] <- 0
      res.coef[q,] <- pmax(pmin(run.coef, ols.max), -ols.max)
    }
    if(SL.method=="nnls"){
      run.coef <- nnls(reg.a, rsp.a.cv[,q])$x
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
  a.SL <- matrix(NA, nrow=dim(a.hat)[1], ncol=Q)
  for(q in 1:Q)
  {
    a.SL[,q] <- a.hat[,q,] %*% matrix(res.coef[q,], ncol=1)
  }
  ###
  return(list(lib.predict=a.hat, SL.predict=a.SL))
}

