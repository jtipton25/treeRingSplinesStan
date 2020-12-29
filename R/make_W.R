#' Generate spline penalty adjacency matrix
#'
#' A function for setting up a conditional autoregressive (CAR) or simultaneous autoregressive (SAR) precision matrix for use as a prior in Bayesian models
#'
#' @param df is the degree of freedom parameter
#' @return A \eqn{df \times df}{df x df} adjacency matrices
#' @importFrom igraph as_adjacency_matrix make_lattice
#' @importFrom Matrix Diagonal colSums
#' @importFrom spam spam
#'
#' @examples
#' df <- 4
#' phi <- 0.8
#' Q <- make_Q(df, phi)
#'
#' @export

make_W <- function(df) {

    if (!is_positive_integer(df, 1))
        stop("df must be a positive integer")

    W <- as_adjacency_matrix(make_lattice(length = df, dim = 2))
    return(as.matrix(W))
}
