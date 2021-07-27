#' Estimating function with augmentation for the parameters of interest in two-phase design.
#' 
#' 'psiStar()' is estimation function with augmentation in the two-phase design.
#' 
#'  This function is used to calculate the estimation function with augmentation in the two-phase design.
#'  
#' @param Y a vector for outcome
#' @param Trt a vector for treatment indicator
#' @param Z a vector for biomarker
#' @param R an indicator for whether the biomarker is observed 
#' @param Pi probability of receiving the treatment (T=1) 
#' @param probR conditional probability for the event that biomarker is observed
#' @param bet a vector containing the parameters of interest
#' @param a.fun the first augmented term which is a function of  (biomarker, baseline covairate)
#' @param b.fun the second augmented term which is a function of (outcome, treatment, baseline covairate)
#'
#'  @return a matrix containing the estimating function at each individual observation, in which
#'  the row represents the subjects, the column represents the parameters
#'  
#'  @export
#'  
psiStar <- function(Y, Trt, Z, R, Pi, probR, bet, a.fun, b.fun){
  
  n <- length(Y)
  datMat <- cbind(Y, Trt, Z)
  Q <- 4
  psi.value <- t(apply(datMat, MARGIN=1, FUN=psi.fun, bet=bet))
  temp1 <- matrix(R/probR, nrow=n, ncol=Q, byrow=F)*psi.value
  temp2 <- matrix(R*(Trt-Pi)/probR, nrow=n, ncol=Q, byrow=F)*a.fun
  temp3 <- matrix((R-probR)/probR, nrow=n, ncol=Q, byrow=F)*b.fun
  res.mat <- temp1-temp2-temp3
  return(res.mat)
}