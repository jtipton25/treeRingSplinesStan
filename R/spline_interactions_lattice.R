#' Generate separable B-spline interactions multivariate regression penalty
#' precision matrix for prior
#'
#' @param df The degree of freedom for the basis function expansion
#'
#' @return A matrix with df^2 rows and 2 columns that indicate which basis
#'     functions are used in the interaction
#' @export
#'
#' @importFrom igraph make_lattice as_adjacency_matrix
#'
#' @examples
#'
#' set.seed(111)
#' df <- 5
#'
#' Q <- spline_interactions_penalty(5)
#'

spline_interactions_lattice <- function(df) {

  if (!is_positive_integer(df, 1))
    stop("df must be a positive integer")



  # lattice <- as_adjacency_matrix(make_lattice(length = df, dim = 2))
  idx <- expand.grid(1:df, 1:df)
  return(idx)
}
