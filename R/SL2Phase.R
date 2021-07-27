#' Efficient Estimation for biomarker studies in the two-phase sampling design.
#' 
#' 'SL2Phase()' returns efficient estimation  for biomarker studies in the two-phase
#'  design in which the biomarker is only measured on a random sub-sample of subjects
#'  chosen on the basis of (Y, T, W).
#'   
#'  This function is used to estimate the parameter for the relationship between a 
#'  clinical outcome of interest and the biomarker, treatment indicator (can be omitted
#'  in the case of a prognostic biomarker) in the two-phase setting. The efficient 
#'  estimation incorporates the independence between treatment and all relevant
#'   baseline covariates. 
#'   
#' @param Y a vector for outcome
#' @param Trt a vector for treatment indicator
#' @param Z a vector for biomarker
#' @param R an indicator for whether the biomarker is observed 
#' @param W a matrix for baseline covariates (not inlcuing the biomarker)
#' @param Pi probability of receiving the treatment (T=1) 
#' @param probR conditional probability for the event that biomarker is observed
#' @param SL.lib a list of functions for candidate prediction algorithms
#' @param SL.family allows gaussian or binomial to describe the error distribution
#' @param SL.method a character value indicating the method used to estimate the coefficients
#' for the super learner and the model to combine the individual algorithms
#' @param ols.max upper bound for the coefficients of individual algorithms in the super learner
#' @param SL.cv  Number of splits for the cross-validation step in super learner
#' @param Split.cv Number of splits for the calculation of cross-validated variance estimate
#' @param iterMax maximum iteration numbers  
#' @param iterTol iteration tolerance
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
SL2Phase <- function(Y, Trt, Z, R, W, Pi, probR, SL.lib, SL.family, SL.method, ols.max, SL.cv, split.cv, iterMax, iterTol){
  
  n <- length(Y)
  Q <- 4
  if(!is.matrix(W))  W <- matrix(W, ncol=1)
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
  
  ######################################################################################
  ##############  Step 1:  Efficient Estimator using iteration
  if(iterMax==0)
  {
    a.hat <- a.Ini
    b.hat <- b.Ini
    bet.hat <- bet.Ini
  }else{
    ##### estimate the functions a and b iteatively 
    a.old <- matrix(0, nrow=n, ncol=Q)
    b.old <- matrix(0, nrow=n, ncol=Q)
    trace.old <- sum(diag(psiStarVar(Y=Y, Trt=Trt, Z=Z, R=R, Pi=Pi, probR=probR, bet=bet.Ini, a.fun=a.old, b.fun=b.old)))
    for(iter in 1:iterMax)
    {
      #### update a() iteratively using the super learner algorithm
      SL.run1 <- SuperLearnerA(Y=Y, Trt=Trt, Z=Z, R=R, W=W, Pi=Pi, probR=probR, newZ=Z, newW=W, SL.lib=SL.lib,
                               SL.family=SL.family, SL.method=SL.method, b.fun=b.old, ols.max=ols.max, SL.cv=SL.cv)
      a.new <- SL.run1$SL.predict
      #### update b() iteratively using the super learner algorithm
      SL.run2 <- SuperLearnerB(Y=Y, Trt=Trt, Z=Z, R=R, W=W, Pi=Pi, probR=probR, newY=Y, newTrt=Trt, newW=W, SL.lib=SL.lib,
                               SL.family=SL.family, SL.method=SL.method, a.fun=a.new, ols.max=ols.max, SL.cv=SL.cv)
      b.new <- SL.run2$SL.predict
      ### caculate the trace difference to determine whether further iteation is needed
      trace.new <- sum(diag(psiStarVar(Y=Y, Trt=Trt, Z=Z, R=R, Pi=Pi, probR=probR, bet=bet.Ini, a.fun=a.new, b.fun=b.new)))
      diff.trace <- trace.old-trace.new
      if(iterMax>2 & abs(diff.trace)<iterTol)
      {
        print(iter) 
        break
      }else{
        ### update a and b
        a.old <- a.new
        b.old <- b.new
        trace.old <- trace.new
      }
    } ### complete the iteration  
    
    
    a.hat <- a.new
    b.hat <- b.new
    #### update the estimate of beta
    psiStar.update <- function(bet)
    {
      return(colSums(psiStar(Y=Y, Trt=Trt, Z=Z, R=R, Pi=Pi, probR=probR, bet=bet, a.fun=a.hat, b.fun=b.hat)))
    }
    bet.hat <- multiroot(f=psiStar.update, start=rep(0, Q), maxiter=500)$root
  }
  
  ######################################################################################
  ########  Step 2: Direct variance estimate for beta #######################
  var.hat <- BetaVar(Y=Y, Trt=Trt, Z=Z, R=R, Pi=Pi, probR=probR, bet=bet.hat, a.fun=a.hat, b.fun=b.hat)
  
  
  ##############################################################################
  ######### Step 3: Cross-validated variance estimate for beta (iterMax > 0)
  ### split the sample randomly into split.cv subsamples that are roughly equal in size
  if(iterMax==0)
  {
    var.hat.SS <- matrix(NA, nrow=Q, ncol=Q)
  }else{
    a.hat.SS <- matrix(NA, nrow=n, ncol=Q)
    b.hat.SS <- matrix(NA, nrow=n, ncol=Q)
    fold.SS <- sample(1:split.cv, n, replace=TRUE)
    for(h in 1:split.cv)
    {
      SS.train <- which(fold.SS!=h)
      SS.valid <- which(fold.SS==h)
      Y.SS <- Y[SS.train]; Trt.SS <- Trt[SS.train]; Z.SS <- Z[SS.train]
      W.SS <- W[SS.train,]; R.SS <- R[SS.train]; probR.SS <- probR[SS.train]
      ### obtain the initial estimate of beta based on the training sample
      aIni.SS <- matrix(0, nrow=length(SS.train), ncol=Q)
      bIni.SS <- matrix(0, nrow=length(SS.train), ncol=Q)
      psiStarIni.SS <- function(bet)
      {
        return(colSums(psiStar(Y=Y.SS, Trt=Trt.SS, Z=Z.SS, R=R.SS, probR=probR.SS, Pi=Pi, bet=bet, a.fun=aIni.SS, b.fun=bIni.SS)))
      }
      betIni.SS <- multiroot(f=psiStarIni.SS, start=rep(0, Q), maxiter=500)$root
      
      ###### iteration
      a.old.SS <- matrix(0, nrow=length(SS.train), ncol=Q)
      b.old.SS <- matrix(0, nrow=length(SS.train), ncol=Q)
      trace.old.SS <- sum(diag(psiStarVar(Y=Y.SS, Trt=Trt.SS, Z=Z.SS, R=R.SS, probR=probR.SS, Pi=Pi, bet=betIni.SS,
                                          a.fun=a.old.SS, b.fun=b.old.SS)))
      for(iter in 1:iterMax)
      {
        #### update a() iteratively using the super learner algorithm
        SL1.SS <- SuperLearnerA(Y=Y.SS, Trt=Trt.SS, Z=Z.SS, R=R.SS, W=W.SS, Pi=Pi, probR=probR.SS, newZ=Z, newW=W, SL.lib=SL.lib,
                                SL.family=SL.family, SL.method=SL.method, b.fun=b.old.SS, ols.max=ols.max, SL.cv=SL.cv)$SL.predict
        a.new.SS <- SL1.SS[SS.train,]
        a.new.valid <- SL1.SS[SS.valid,]
        #### update b() iteratively using the super learner algorithm
        SL2.SS <- SuperLearnerB(Y=Y.SS, Trt=Trt.SS, Z=Z.SS, R=R.SS, W=W.SS, Pi=Pi, probR=probR.SS, newY=Y, newTrt=Trt, newW=W, SL.lib=SL.lib, 
                                SL.family=SL.family, SL.method=SL.method, a.fun=a.new.SS, ols.max=ols.max, SL.cv=SL.cv)$SL.predict
        b.new.SS <- SL2.SS[SS.train,]
        b.new.valid <- SL2.SS[SS.valid,]
        ### caculate the trace difference to determine whether further iteation is needed
        trace.new.SS <- sum(diag(psiStarVar(Y=Y.SS, Trt=Trt.SS, Z=Z.SS, R=R.SS, probR=probR.SS, Pi=Pi, bet=betIni.SS, 
                                            a.fun=a.new.SS, b.fun=b.new.SS)))
        diff.trace.SS <-  trace.old.SS-trace.new.SS
        if(iterMax > 2 & abs(diff.trace.SS) < iterTol)
        {
          print(iter)  
          break
        }else{
          ### update a and b
          a.old.SS <- a.new.SS
          b.old.SS <- b.new.SS
          trace.old.SS <- trace.new.SS
        }
      }
      a.hat.SS[SS.valid,] <- a.new.valid
      b.hat.SS[SS.valid,] <- b.new.valid
    }
    #### update the estimate of beta
    psiStar.SS <- function(bet)
    {
      return(colSums(psiStar(Y=Y, Trt=Trt, Z=Z, R=R, Pi=Pi, probR=probR, bet=bet, a.fun=a.hat.SS, b.fun=b.hat.SS)))
    }
    bet.hat.SS <- multiroot(f=psiStar.SS, start=rep(0, Q), maxiter=500)$root
    var.hat.SS <- BetaVar(Y=Y, Trt=Trt, Z=Z, R=R, Pi=Pi, probR=probR, bet=bet.hat.SS, a.fun=a.hat.SS, b.fun=b.hat.SS)
  }
  
  ################################################################
  ### return the results (point estimate and variance estimate of beta)
  return(list(pEst=bet.hat, varEst.cv=var.hat.SS))
}

