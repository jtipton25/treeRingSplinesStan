#' Generate separable B-spline interactions with multivariate regression penalty
#'
#' @param x1 A vector of length n representing the first variable to use to make the basis function expansion
#' @param x2 A vector of length n representing the second variable to use to make the basis function expansion
#' @param df The degree of freedom for the basis function expansion
#'
#' @return A matrix of separable B-spline interactions
#' @export
#'
#' @importFrom igraph make_lattice
#' @importFrom splines bs
#'
#' @examples
#'
#' library(splines)
#' set.seed(111)
#' X <- rnorm(10)
#' Y <- rnorm(10)
#'
#' Z1 <- spline_interactions(X, Y, 5)
#' Z2 <- model.matrix( ~ bs(X, df = 5, intercept = FALSE):bs(Y, df = 5, intercept = FALSE) - 1)
#'
#' all.equal(Z1, Z2, check.attributes = FALSE)
#'
spline_interactions <- function(x1, x2, df) {


  ## if making into a package, use my predicates for testing input
  ## also make unit tests
  if (!is_numeric_vector(x1, length(x1)))
    stop("x1 must be a numeric vector")
  if (!is_numeric_vector(x2, length(x2)))
    stop("x2 must be a numeric vector")
  if (length(x1) != length(x2))
    stop("x1 and x2 must each be of the same length")
  if (!is_positive_integer(df, 1))
    stop("df must be a positive integer")

  n <- length(x1)

  X <- matrix(0, n, df^2)

  x1_bs <- bs(x1, df = df, intercept = FALSE)
  x2_bs <- bs(x2, df = df, intercept = FALSE)

  idx <- expand.grid(1:df, 1:df)

  ## could make this more efficient -- check example on desktop
  for (i in 1:(df^2)) {
    X[, i] <-  x1_bs[, idx[i, 1]] * x2_bs[, idx[i, 2]]
  }

  return(X)
}


