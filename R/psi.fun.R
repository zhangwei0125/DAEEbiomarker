#' @title Estimating function without augmentation for the parameters of interest.
#'
#' @description 'psi.fun()' is the value of the estimating function for the individual observation.
#'
#' @details This function is used to calculate the estimating function for the individual observation.
#'
#' @param dat.vec a vector contianing the observation of outcome, treatment, biomaker
#' @param bet value for the parameter beta
#'
#' @return a value of the estimating function for the individual observation (outcome, treatment, biomaker)
#'
#' @example
#'  dat.vec = c(1.5, 1, 0.3)
#'  bet = 0.5
#'  psi.fun(dat.vec=dat.vec, bet=bet)
#'
#' @export
#'
psi.fun <- function(dat.vec, bet){

  y <- dat.vec[1] ##outcome
  trt <- dat.vec[2] ##treatment
  z <- dat.vec[-c(1,2)] ##biomarker
  comVec <- c(1, trt, z, trt*z)
  linearSum <- c(matrix(comVec, nrow=1)%*%matrix(bet, ncol=1))
  res.vec <- (y-invLogit(linearSum))*matrix(comVec, ncol=1)
  return(res.vec)
}
