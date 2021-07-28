#' @title The first derivative of the inverse of the logistic function.
#'
#' @description 'invLogitD1()' returns the first derivative value of the inverse of the logistic function.
#'
#' @details This function is used to calculate the first derivative function for the inverse logit function.
#'
#' @param u a real value
#'
#' @return  the first derivative value
#'
#' @examples
#' invLogitD1(2)
#'
#' @export
#'
invLogitD1 <- function(u){

  return(exp(u)/(1+exp(u))^2)
}
