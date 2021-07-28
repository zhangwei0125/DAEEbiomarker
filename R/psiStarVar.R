#' @title Variance estimate for the estimating function with augmentation in the two-phase design.
#'
#' @description 'psiStarVar()' returns the variance estimate of the estimating function with augmentation in the two-phase design.
#'
#' @details This function is used to calculate the variance of the augmented estimating function
#' in the two-phase design.
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
#'  @return a matrix for the variance estimate of the estimating function with augmentation
#'  in the two-phase deisgn
#'
#'  @export
#'
psiStarVar <- function(Y, Trt, Z, R, Pi, probR, bet, a.fun, b.fun){

  n <- length(Y)
  psi.star <- psiStar(Y=Y, Trt=Trt, Z=Z, R=R, Pi=Pi, probR=probR, bet=bet, a.fun=a.fun, b.fun=b.fun)
  res.var <- matrix(0, nrow=ncol(psi.star), ncol=ncol(psi.star))
  for(i in 1:n)
  {
    res.var <- res.var + 1/n*(matrix(psi.star[i,], ncol=1) %*% matrix(psi.star[i,], nrow=1))
  }
  ################################################################
  return(res.var)
}
