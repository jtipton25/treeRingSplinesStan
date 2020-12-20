#' Generate separable B-spline interactions with multivariate regression penalty
#'
#' @param X A $n x p$ matrix of covariates from which you want to evaluate
#'    all pairwise combinations of b-spline basis interactions.
#' @param df The degree of freedom for the basis function expansion
#'
#' @return An array of separable B-spline interactions. The first dimension is
#'    the number of pairwise interactions (p choose 2). The second dimension is the number
#'    of observations $n$ and the third dimension is the number of knots in the
#'    basis expansion $df^2$.
#' @export
#'
#' @importFrom igraph make_lattice
#' @importFrom splines bs
#' @importFrom utils combn
#'
#' @examples
#'
#' library(splines)
#' set.seed(111)
#' X <- matrix(rnorm(40), 10, 4) ## 10 observations and 4 variables
#'
#' Z <- build_spline_interactions(X, 5)
#'
build_spline_interactions <- function(X, df) {


  ## if making into a package, use my predicates for testing input
  ## also make unit tests
  if (!is_numeric_matrix(X, nrow(X), ncol(X)))
    stop("X must be a numeric matrix with at least 2 columns")
  if (ncol(X) < 2)
    stop("X must be a numeric matrix with at least 2 columns")

  if (!is_positive_integer(df, 1))
    stop("df must be a positive integer")

  n <- nrow(X)
  p <- ncol(X)

  if (p > 20)
    stop("Currently, we only accommodate p <= 20 covariate values")

  n_interactions <- choose(p, 2)

  ## check if columnnames are specified
  if (is.null(colnames(X)))
    colnames(X) <- paste("Var", 1:p, sep = "-")

  interaction_names = apply(t(combn(colnames(X), 2)), 1, paste, collapse = ":")

  Z <- array(0,
             dim = c(n_interactions, n, df^2),
             dimnames = list(
               variables   = interaction_names,
               observation = 1:n,
               knot        = 1:df^2
             )
  )


  interactions <- t(combn(p, 2))
  ## mapping the product of basis functions
  idx <- expand.grid(1:df, 1:df)

  ## pre-compute the basis functions
  X_bs <- array(0, dim = c(p, n, df))
  for (i in 1:p) {
    X_bs[i, , ] <- bs(X[, i], df = df, intercept = FALSE)
  }

  for (i in 1:n_interactions) {
    for (j in 1:(df^2)) {
      Z[i, , j] <- X_bs[interactions[i, 1], , idx[j, 1]] * X_bs[interactions[i, 2], , idx[j, 2]]
    }
  }

  # idx <- expand.grid(1:df, 1:df)
  #
  # ## could make this more efficient -- check example on desktop
  # for (i in 1:(df^2)) {
  #   X[, i] <-  x1_bs[, idx[i, 1]] * x2_bs[, idx[i, 2]]
  # }

  return(Z)
}


