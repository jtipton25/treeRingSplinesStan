
#' lm_linear_no_p
#'
#' @param y A vector of length n of tree ring increments on a log scale
#' @param X A n times p matrix of covariate values
#' @param X_int A n times n_tree matrix of covariate values
#' @param X_pred A n_pred by p matrix of covariate values at which to generate predictions
#' @param ... Additional arguments to `rstan::sampling()`
#'
#' @return A stan object that is the result of fitting the linear regression model
#' @export
#'
lm_linear_X <- function(y, X, X_int, X_pred, X_int_pred, ...) {
  ## add in error checking and unit testing

  standata <- list(
    y          = y,
    X          = X,
    X_int      = X_int,
    n          = length(y),
    K          = ncol(X),
    n_tree     = n_tree,
    n_pred     = nrow(X_pred),
    X_pred     = X_pred,
    X_int_pred = X_int_pred
  )

  out <- rstan::sampling(stanmodels$linear_X, data = standata, ...)
  return(out)

}


