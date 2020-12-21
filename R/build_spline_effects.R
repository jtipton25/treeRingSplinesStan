#' Build the non-linear response from the basis functions and paramters
#'
#' @param Z A q by n by df array of spline-expanded coefficients
#' @param beta A q by df matrix of spline coefficieints
#'
#' @return A n by q matrix of effects (non-linear functional response)
#' @export
#'
build_spline_effects <- function(Z, beta) {

  # check that the array is numeric
  # if(!is_numeric_array(Z, dim(Z)))
  #   stop("Z must by a numeric array with three dimesions")
  # if(length(dim(Z)) != 3)
  #   stop("Z must by a numeric array with three dimesions")
  # check the dimensions
  if (dim(Z)[1] != dim(beta)[1])
    stop("The number of rows of Z must be equal to the number of rows of beta")
  if (dim(Z)[3] != dim(beta)[2])
    stop("The number of slices of Z (3rd dimension) must be equal to the number of columns of beta")

  q <- dim(Z)[1]
  n <- dim(Z)[2]
  effects <- sapply(1:q, function(i) Z[i, , ] %*% beta[i, ])

  return(effects)

}
