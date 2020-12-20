#' Generate separable B-spline interactions multivariate regression penalty
#' precision matrix for prior
#'
#' @param df The degree of freedom for the basis function expansion
#' @param phi is a number between -1 and 1 that defines the strength of the autoregressive process. Typically this will be set to 1 for use as a prior in penalized Bayesian models
#'
#' @return A precision matrix for the smoothing penalty prior
#' @export
#'
#' @examples
#'
#' set.seed(111)
#' df <- 5
#'
#' Q <- spline_interactions_penalty(5)
#'
spline_interactions_penalty <- function(df, phi = 1) {

  if (!is_positive_integer(df, 1))
    stop("df must be a positive integer")

  if (!is_numeric(phi, 1))
    stop("phi must be a number with values between -1 and 1.")
  if (phi < -1 | phi > 1)
    stop("phi must be a number with values between -1 and 1.")

  Q <- as.matrix(make_Q(df, phi = phi))

  return(Q)
}


