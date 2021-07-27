#' First derivative of the estimating function without augmentation at a certain value. 
#' 
#' 'Der1.psi()' is the first derivative value of the estimating function for the
#'  individual observation at a certain value. 
#' 
#'  This function is used to calculate the first derivative of the estimating function
#'   for the individual observation at a certain value.
#'  
#'  @param dat.vec a vector contianing the observation of outcome, treatment, biomaker
#'  @param bet value for the parameter beta
#'  
#'  @return a value of the first derivative of the estimating function for the individual
#'   observation (outcome, treatment, biomaker)
#'  
#'  @example 
#'  dat.vec = c(1.5, 1, 0.3)
#'  bet = 0.5
#'  Der1.psi(dat.vec=dat.vec, bet=bet)
#'  
#'  @export
Der1.psi <- function(dat.vec, bet){
  
  y <- dat.vec[1]
  trt <- dat.vec[2]
  z <- dat.vec[-c(1,2)]
  comVec <- c(1, trt, z, trt*z)
  linearSum <- c(matrix(comVec, nrow=1)%*%matrix(bet, ncol=1))
  res.mat <-  (-invLogitD1(linearSum))*(matrix(comVec, ncol=1) %*% matrix(comVec, nrow=1))
  return(res.mat)
}