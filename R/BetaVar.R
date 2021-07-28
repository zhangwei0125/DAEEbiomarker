#' @title Sandwich variance estimation for the biomarker studies in the two-phase setting.
#'
#' @description 'BetaVar()' returns the variance estimate of the parameter beta in the two-phase setting.
#'
#' @details This function is used to calculate the sandwich variance estimate for the parameter estimate
#' in the two-phase setting.
#'
#' @param Y a vector for outcome
#' @param Trt a vector for treatment indicator
#' @param Z a vector for biomarker
#' @param R an indicator for whether the biomarker is observed
#' @param Pi probability of receiving the treatment (T=1)
#' @param probR conditional probability for the event that biomarker is observed
#' @param bet parameter estimate
#' @param a.fun the first augmented term which is a function of (biomarker, baseline covairate)
#' @param b.fun the first augmented term which is a function of (outcome, treatment, baseline covairate)
#'
#' @return the variance estimate of the parameter beta
#'
#' @export
#'
BetaVar <- function(Y, Trt, Z, R, Pi, probR, bet, a.fun, b.fun){

  n <- length(Y)
  datMat <- cbind(Y, Trt, Z)
  Q <- 4
  ### calculate the psi.star
  psi.star <- psiStar(Y=Y, Trt=Trt, Z=Z, R=R, Pi=Pi, probR=probR, bet=bet, a.fun=a.fun, b.fun=b.fun)
  temp.D <- matrix(0, Q, Q)
  temp.v <- matrix(0, Q, Q)
  for(i in 1:n)
  {
    temp.D <- temp.D + 1/n*(R[i]/probR[i])*Der1.psi(dat.vec=datMat[i,], bet=bet)
    temp.v <- temp.v + 1/n*(matrix(psi.star[i,], ncol=1)%*%matrix(psi.star[i,], nrow=1))
  }
  temp.invD <- solve(temp.D)
  res.var <- 1/n*(temp.invD%*%temp.v%*%t(temp.invD))
  ################################################################
  return(res.var)
}


