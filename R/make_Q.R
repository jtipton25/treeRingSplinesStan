#' Generate spline penalty precision matrix
#'
#' A function for setting up a conditional autoregressive (CAR) or simultaneous autoregressive (SAR) precision matrix for use as a prior in Bayesian models
#'
#' @param df is the degree of freedom parameter
#' @param phi is a number between -1 and 1 that defines the strength of the autoregressive process. Typically this will be set to 1 for use as a prior in penalized Bayesian models
#' @param prec_model is a string that takes the values "CAR" or "SAR" and defines the graphical structure for the precision matrix.
#' @return a list of  \eqn{n \times n}{n x n} sparse spam matrices or Matrix matrices of class "dgCMatrix" (see Matrix package for details)
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

make_Q <- function(df, phi, prec_model = "CAR") {

    if (!is_numeric(phi, 1))
        stop("phi must be a number with values between -1 and 1.")
    if (phi < -1 | phi > 1)
        stop("phi must be a number with values between -1 and 1.")

    if (!is_positive_integer(df, 1))
        stop("df must be a positive integer")

    if (!(prec_model %in% c("CAR", "SAR")))
        stop('The only valid options for prec_model are "CAR" and "SAR".')

    Q <- NULL
    W <- as_adjacency_matrix(make_lattice(length = df, dim = 2))
    D <- Diagonal(x = colSums(W))
    if (prec_model == "CAR") {
        Q <- D - phi * W
    } else if (prec_model == "SAR") {
        B <- diag(nrow(W)) - phi * W %*% Diagonal(x = 1 / colSums(W))
        Q <- t(B) %*% B
    }
    Q <- as.matrix(Q)

    return(Q)
}
