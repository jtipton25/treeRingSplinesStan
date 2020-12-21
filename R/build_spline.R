#' Generate univariate B-spline expansions
#'
#' @param X A $n x p$ matrix of covariates from which you want to evaluate
#'    all marginal spline effects.
#' @param df The degree of freedom for the basis function expansion
#'
#' @return An array of univariate B-spline function expansions. The first dimension is
#'    the number of variables $p$. The second dimension is the number
#'    of observations $n$ and the third dimension is the number of knots in the
#'    basis expansion $df$.
#' @export
#'
#' @importFrom splines bs
#'
#' @examples
#'
#' library(splines)
#' set.seed(111)
#' X <- matrix(rnorm(40), 10, 4) ## 10 observations and 4 variables
#'
#' Z <- build_spline(X, 5)
#'
build_spline <- function(X, df) {


  ## if making into a package, use my predicates for testing input
  ## also make unit tests
  if (!is_numeric_matrix(X, nrow(X), ncol(X)))
    stop("X must be a numeric matrix with at least 1 columns")
  if (ncol(X) < 1)
    stop("X must be a numeric matrix with at least 1 columns")

  if (!is_positive_integer(df, 1))
    stop("df must be a positive integer")

  n <- nrow(X)
  p <- ncol(X)

  if (p > 20)
    stop("Currently, we only accommodate p <= 20 covariate values")

  ## check if colnames are specified
  if (is.null(colnames(X)))
    colnames(X) <- paste("Var", 1:p, sep = "-")


  ## pre-compute the basis functions
  Z <- array(0, dim = c(p, n, df))
  for (i in 1:p) {
    Z[i, , ] <- bs(X[, i], df = df, intercept = FALSE)
  }


  return(Z)
}


