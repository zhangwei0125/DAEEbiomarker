#' The inverse of the logistic function.
#' 
#' 'invLogit()' returns a value in the interval [0,1].
#' 
#' This function is used to transform a value in the real line to the interval [0,1].
#' 
#' @param u a real value
#' 
#' @return a value in the interval [0,1]
#' 
#' @examples
#' invLogit(2)
#' 
#' @export
invLogit <- function(u){
  
  return(exp(u)/(1+exp(u)))
}