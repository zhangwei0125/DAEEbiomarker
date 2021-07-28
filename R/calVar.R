#' @title Sandwich variance estimation for the biomarker studies in the complete-data setting.
#'
#' @description 'calVar()' returns the variance estimate of the parameter beta in the complete-data setting.
#'
#' @details This function is used to calculate the sandwich variance estimate for the parameter estimate
#' in the complete-data setting.
#'
#' @param Y a vector for outcome
#' @param Trt a vector for treatment indicator
#' @param Z a vector for biomarker
#' @param Pi probability of receiving the treatment (T=1)
#' @param bet parameter estimate
#' @param a.fun augmented function evaluated at the observation (Y, Trt, Z)
#'
#' @return the variance estimate of the parameter beta
#'
#' @export
#'
calVar <- function(Y, Trt, Z, Pi, bet, a.fun){

  datMat <- cbind(Y, Trt, Z)
  Q <- 4
  n <- length(Y)
  temp.D <- matrix(0, Q, Q)
  temp.v <- matrix(0, Q, Q)
  for(i in 1:n)
  {
    temp.D <- temp.D + 1/n*Der1.psi(dat.vec=datMat[i,], bet=bet)
    temp.psi <- psi.fun(dat.vec=datMat[i,], bet=bet)-(Trt[i]-Pi)*a.fun[i,]
    temp.v <- temp.v + 1/n*(matrix(temp.psi, ncol=1)%*%matrix(temp.psi, nrow=1))
  }
  temp.invD <- solve(temp.D)
  var.est <- (temp.invD%*%temp.v%*%t(temp.invD))/n
  ################################################################
  return(var.est)
}

